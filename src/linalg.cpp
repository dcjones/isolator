
#include "common.hpp"
#include "config.h"
#include "cpuid.hpp"
#include "linalg.hpp"
#include "logger.hpp"

#include <math.h>
#include <immintrin.h>


/* The implementation of each of these functions is determined at run time */
void* (*aalloc)(size_t n) = NULL;
void (*afree)(void*) = NULL;
void (*acopy)(void* dest, const void* src, size_t n) = NULL;
float (*dotlog)(const float* xs, const float* ys, const size_t n) = NULL;
float (*dotlogc)(const float* xs, const float* ys, const size_t n, const float c) = NULL;
void (*asxpy)(float* xs, const float* ys, const float c,
              const unsigned int* idx, const unsigned int off,
              const size_t n) = NULL;
float (*asxtydsz)(const float* xs, const float* ys, const float* zs,
                  const unsigned int* idx, const unsigned int off,
                  const size_t n) = NULL;
float (*dot)(const float* xs, const float* ys, const float* zs, size_t n);
static __m128 (*log2_sse)(__m128 x) = NULL;
const char* LINALG_INSTR_SET = "";


#ifdef _MSC_VER
  #define ALIGN32_START __declspec(align(32))
  #define ALIGN32_END
  #define ALIGN16_START __declspec(align(16))
  #define ALIGN16_END
#else
  #define ALIGN32_START
  #define ALIGN32_END __attribute__((aligned(32)))
  #define ALIGN16_START
  #define ALIGN16_END __attribute__((aligned(16)))
#endif

#define PS32_CONST(name, c) \
    static const ALIGN32_START float ps32_##name[8] ALIGN32_END = {c, c, c, c, c, c, c, c}

#define PI32_CONST(name, c) \
    static const ALIGN32_START int pi32_##name[8] ALIGN32_END = {c, c, c, c, c, c, c, c}

#define PS16_CONST(name, c) \
    static const ALIGN16_START float ps16_##name[8] ALIGN16_END = {c, c, c, c}

#define PI16_CONST(name, c) \
    static const ALIGN16_START int pi16_##name[8] ALIGN16_END = {c, c, c, c}


float fastlog2(float x_)
{
    if (x_ <= 0.0f) return -INFINITY;

    union {
        float   f;
        int32_t i;
    } x, y, one;

    one.f = 1.0f;
    x.f = x_;

    float e = (float) ((x.i >> 23) - 127);

    y.i = (x.i & 0x007FFFFF) | one.i;
    float m = y.f;

    float p = -3.4436006e-2f;
    p = m * p +  3.1821337e-1f;
    p = m * p + -1.2315303f;
    p = m * p +  2.5988452f;
    p = m * p + -3.3241990f;
    p = m * p +  3.1157899f;

    p *= m - 1.0;

    return p + e;
}


float fastlog(float x)
{
    return fastlog2(x) / M_LOG2E;
}


/* AVX versions */

PS32_CONST(1, 1.0f);

//union {
    //uint32_t i;
    //float    f;
//} inv_mant_mask_int = {~0x7f800000};
//const unsigned int inv_mant_mask_int = ~0x7f800000;
//PI32_CONST(inv_mant_mask, *reinterpret_cast<const float*>(&inv_mant_mask_int));
PI32_CONST(inv_mant_mask, ~0x7f800000);
PI16_CONST(0x7f, 0x7f);

PS32_CONST(log2_c0, 3.1157899f);
PS32_CONST(log2_c1, -3.3241990f);
PS32_CONST(log2_c2, 2.5988452f);
PS32_CONST(log2_c3, -1.2315303f);
PS32_CONST(log2_c4, 3.1821337e-1f);
PS32_CONST(log2_c5, -3.4436006e-2f);
PS32_CONST(neginf, (float) -INFINITY);


void* aalloc_avx(size_t n)
{
    void* xs = _mm_malloc(n, 32);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree_avx(void* xs)
{
    _mm_free(xs);
}


void acopy_avx(void* dest_, const void* src_, size_t n)
{
    assert(n % 4 == 0);
    n /= 4; /* bytes to 4-byte words */

    float* dest = reinterpret_cast<float*>(dest_);
    const float* src = reinterpret_cast<const float*>(src_);

    __m256 x;
    size_t i;
    for (i = 0; i < n / 8; ++i) {
        x = _mm256_load_ps(src + 8 * i);
        _mm256_store_ps(dest + 8 * i, x);
    }

    i *= 8;
    switch (n % 8) {
        case 7: dest[i] = src[i]; ++i;
        case 6: dest[i] = src[i]; ++i;
        case 5: dest[i] = src[i]; ++i;
        case 4: dest[i] = src[i]; ++i;
        case 3: dest[i] = src[i]; ++i;
        case 2: dest[i] = src[i]; ++i;
        case 1: dest[i] = src[i];
    }
}


/* Hopefully a cleaner version of avx_log2, following intels' amaths. */
static __m256 log2_avx(__m256 x)
{
    /* extract the exponent */
    union {
        __m256i a;
        __m128i b[2];
    } xi;
    xi.a = _mm256_castps_si256(x);
    xi.b[0] = _mm_sub_epi32(_mm_srli_epi32(xi.b[0], 23), *(__m128i*) pi16_0x7f);
    xi.b[1] = _mm_sub_epi32(_mm_srli_epi32(xi.b[1], 23), *(__m128i*) pi16_0x7f);
    __m256 e = _mm256_cvtepi32_ps(xi.a);

    /* extract the mantissa */
    __m256 m = _mm256_and_ps(x, *(__m256*) pi32_inv_mant_mask);
    m = _mm256_or_ps(m, *(__m256*) ps32_1);

    __m256 a = *(__m256*) ps32_log2_c5;
    a = _mm256_add_ps(_mm256_mul_ps(m, a), *(__m256*) ps32_log2_c4);
    a = _mm256_add_ps(_mm256_mul_ps(m, a), *(__m256*) ps32_log2_c3);
    a = _mm256_add_ps(_mm256_mul_ps(m, a), *(__m256*) ps32_log2_c2);
    a = _mm256_add_ps(_mm256_mul_ps(m, a), *(__m256*) ps32_log2_c1);
    a = _mm256_add_ps(_mm256_mul_ps(m, a), *(__m256*) ps32_log2_c0);

    a = _mm256_add_ps(_mm256_mul_ps(a, _mm256_sub_ps(m, *(__m256*) ps32_1)), e);

    /* Handle the log(x) = -INFINITY, for x <= 0 cases. */
    a = _mm256_blendv_ps(a, *(__m256*) ps32_neginf,
                         _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LE_OQ));
    return a;
}


float dotlog_avx(const float* xs, const float* ys, const size_t n)
{
    union {
        __m256 v;
        float f[8];
    } ans;
    ans.v = _mm256_setzero_ps();

    __m256 xv, yv;
    size_t i;
    for (i = 0; i < n / 8; ++i) {
        xv = _mm256_load_ps(xs + 8 * i);
        yv = _mm256_load_ps(ys + 8 * i);
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(xv, log2_avx(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3] +
                 ans.f[4] + ans.f[5] + ans.f[6] + ans.f[7];

    /* handle any overhang */
    i *= 8;
    switch (n % 8) {
        case 7: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 6: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 5: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 4: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 3: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 2: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 1: fans += xs[i] * fastlog2(ys[i]);
    }

    return fans;
}


float dotlogc_avx(const float* xs, const float* ys, const size_t n, const float c)
{
    union {
        __m256 v;
        float f[8];
    } ans;
    ans.v = _mm256_setzero_ps();

    __m256 cv = _mm256_set1_ps(c);

    __m256 xv, yv;
    size_t i;
    for (i = 0; i < n / 8; ++i) {
        xv = _mm256_load_ps(xs + 8 * i);
        yv = _mm256_mul_ps(cv, _mm256_load_ps(ys + 8 * i));
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(xv, log2_avx(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3] +
                 ans.f[4] + ans.f[5] + ans.f[6] + ans.f[7];

    /* handle any overhang */
    i *= 8;
    switch (n % 8) {
        case 7: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 6: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 5: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 4: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 3: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 2: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 1: fans += xs[i] * fastlog2(c * ys[i]);
    }

    return fans;
}


void asxpy_avx(float* xs, const float* ys, const float c,
               const unsigned int* idx, const unsigned int off,
               const size_t n)
{
    __m256 yv;
    __m256 cv = _mm256_set1_ps(c);
    union {
        __m256 v;
        float f[8];
    } x;

    union {
        __m128i a;
        unsigned int b[4];
    } i0, i1;

    __m128i voff = _mm_set1_epi32((int) off);

    size_t i;
    for (i = 0; i < 8 * (n / 8); i += 8) {
        yv = _mm256_mul_ps(cv, _mm256_load_ps(ys + i));

        i0.a = _mm_load_si128((__m128i*)(idx + i));
        i0.a = _mm_sub_epi32(i0.a, voff);

        i1.a = _mm_load_si128((__m128i*)(idx + i + 4));
        i1.a = _mm_sub_epi32(i1.a, voff);

        /* load from xs */
        x.f[0] = xs[i0.b[0]];
        x.f[1] = xs[i0.b[1]];
        x.f[2] = xs[i0.b[2]];
        x.f[3] = xs[i0.b[3]];
        x.f[4] = xs[i1.b[0]];
        x.f[5] = xs[i1.b[1]];
        x.f[6] = xs[i1.b[2]];
        x.f[7] = xs[i1.b[3]];

        x.v = _mm256_add_ps(x.v, yv);

        /* store in xs */
        xs[i0.b[0]] = x.f[0];
        xs[i0.b[1]] = x.f[1];
        xs[i0.b[2]] = x.f[2];
        xs[i0.b[3]] = x.f[3];
        xs[i1.b[0]] = x.f[4];
        xs[i1.b[1]] = x.f[5];
        xs[i1.b[2]] = x.f[6];
        xs[i1.b[3]] = x.f[7];
    }

    /* handle overhang */
    i = 8 * (n / 8);
    switch (n % 8) {
        case 7: xs[idx[i] - off] += c * ys[i]; ++i;
        case 6: xs[idx[i] - off] += c * ys[i]; ++i;
        case 5: xs[idx[i] - off] += c * ys[i]; ++i;
        case 4: xs[idx[i] - off] += c * ys[i]; ++i;
        case 3: xs[idx[i] - off] += c * ys[i]; ++i;
        case 2: xs[idx[i] - off] += c * ys[i]; ++i;
        case 1: xs[idx[i] - off] += c * ys[i];
    }
}


float asxtydsz_avx(const float* xs, const float* ys, const float* zs,
                   const unsigned int* idx, const unsigned int off,
                   const size_t n)
{
    union {
        __m256 v;
        float f[8];
    } ans, x, z;
    ans.v = _mm256_setzero_ps();

    union {
        __m128i a;
        unsigned int b[4];
    } i0, i1;

    __m128i voff = _mm_set1_epi32((int) off);

    size_t i;
    for (i = 0; i < 8 * (n / 8); i += 8) {
        __m256 y = _mm256_load_ps(ys + i);

        i0.a = _mm_load_si128((__m128i*)(idx + i));
        i0.a = _mm_sub_epi32(i0.a, voff);

        i1.a = _mm_load_si128((__m128i*)(idx + i + 4));
        i1.a = _mm_sub_epi32(i1.a, voff);

        /* load from xs */
        x.f[0] = xs[i0.b[0]];
        x.f[1] = xs[i0.b[1]];
        x.f[2] = xs[i0.b[2]];
        x.f[3] = xs[i0.b[3]];
        x.f[4] = xs[i1.b[0]];
        x.f[5] = xs[i1.b[1]];
        x.f[6] = xs[i1.b[2]];
        x.f[7] = xs[i1.b[3]];

        /* load from zs */
        z.f[0] = zs[i0.b[0]];
        z.f[1] = zs[i0.b[1]];
        z.f[2] = zs[i0.b[2]];
        z.f[3] = zs[i0.b[3]];
        z.f[4] = zs[i1.b[0]];
        z.f[5] = zs[i1.b[1]];
        z.f[6] = zs[i1.b[2]];
        z.f[7] = zs[i1.b[3]];

        ans.v = _mm256_add_ps(ans.v, _mm256_div_ps(_mm256_mul_ps(x.v, y), z.v));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3] +
                 ans.f[4] + ans.f[5] + ans.f[6] + ans.f[7];


    /* handle overhang */
    i = 8 * (n / 8);
    switch (n % 8) {
        case 7: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 6: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 5: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 4: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 3: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 2: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 1: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off];
    }

    return fans;
}


float dot_avx(const float* xs, const float* ys, const float* zs, size_t n)
{
    union ans_t
    {
        __m256  v;
        float   f[8];
    } ans;
    ans.v = _mm256_setzero_ps();

    size_t i;
    __m256 x, y, z;
    for (i = 0; i < 8 * (n / 8); i += 8) {
        x = _mm256_load_ps(xs + i);
        y = _mm256_loadu_ps(ys + i);
        z = _mm256_load_ps(zs + i);
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(_mm256_mul_ps(x, y), z));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3] +
                 ans.f[4] + ans.f[5] + ans.f[6] + ans.f[7];

    /* handle overhang */
    i = 8 * (n / 8);
    switch (n % 8) {
        case 7: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 6: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 5: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 4: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 3: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 2: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 1: fans += xs[i] * ys[i] * zs[i];
    }

    return fans;
}


void* aalloc_sse(size_t n)
{
    void* xs = _mm_malloc(n, 16);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree_sse(void* xs)
{
    _mm_free(xs);
}


void acopy_sse(void* dest_, const void* src_, size_t n)
{
    assert(n % 4 == 0);
    n /= 4; /* bytes to 4-byte words */

    float* dest = reinterpret_cast<float*>(dest_);
    const float* src = reinterpret_cast<const float*>(src_);

    __m128 x;
    size_t i;
    for (i = 0; i < n / 4; ++i) {
        x = _mm_load_ps(src + 4 * i);
        _mm_store_ps(dest + 4 * i, x);
    }

    i *= 4;
    switch (n % 4) {
        case 3: dest[i] = src[i]; ++i;
        case 2: dest[i] = src[i]; ++i;
        case 1: dest[i] = src[i];
    }
}


PS16_CONST(1, 1.0f);
PI16_CONST(inv_mant_mask, ~0x7f800000);

PS16_CONST(log2_c0, 3.1157899f);
PS16_CONST(log2_c1, -3.3241990f);
PS16_CONST(log2_c2, 2.5988452f);
PS16_CONST(log2_c3, -1.2315303f);
PS16_CONST(log2_c4, 3.1821337e-1f);
PS16_CONST(log2_c5, -3.4436006e-2f);
PS16_CONST(neginf, (float) -INFINITY);


/* Comptue log2 over an sse single vector. */
static __m128 log2_sse2(__m128 x)
{
    /* extract exponent */
    const __m128i i = _mm_castps_si128(x);
    __m128 e = _mm_cvtepi32_ps(
                 _mm_sub_epi32(
                   _mm_srli_epi32(i, 23),
                   _mm_set1_epi32(127)));

    /* extract mantissa */
    const __m128i mant_mask = _mm_set1_epi32(0x007FFFFF);
    const __m128 one = _mm_set1_ps(1.0f);
    __m128 m = _mm_or_ps(_mm_castsi128_ps(_mm_and_si128(i, mant_mask)), one);

    /* polynomial approximation on the mantissa */
    __m128 p = *(__m128*) ps16_log2_c5;
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c4);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c3);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c2);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c1);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c0);

    p = _mm_add_ps(_mm_mul_ps(p, _mm_sub_ps(m, *(__m128*) ps16_1)), e);

    /* handle the case with x= < 0.0 */
    // TODO: test this
    __m128 mask = _mm_cmpgt_ps(x, _mm_setzero_ps());
    p = _mm_or_ps(_mm_and_ps(mask, p), _mm_andnot_ps(mask, _mm_set1_ps(-INFINITY)));

    return p;
}


static __m128 log2_sse4(__m128 x)
{
    /* extract exponent */
    const __m128i i = _mm_castps_si128(x);
    __m128 e = _mm_cvtepi32_ps(
                 _mm_sub_epi32(
                   _mm_srli_epi32(i, 23),
                   _mm_set1_epi32(127)));

    /* extract mantissa */
    __m128 m = _mm_and_ps(x, *(__m128*) pi16_inv_mant_mask);
    m = _mm_or_ps(m, *(__m128*) ps16_1);

    /* polynomial approximation on the mantissa */
    __m128 p = *(__m128*) ps16_log2_c5;
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c4);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c3);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c2);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c1);
    p = _mm_add_ps(_mm_mul_ps(m, p), *(__m128*) ps16_log2_c0);

    p = _mm_add_ps(_mm_mul_ps(p, _mm_sub_ps(m, *(__m128*) ps16_1)), e);

    /* handle the case with x= < 0.0 */
    p = _mm_blendv_ps(p, _mm_set1_ps(-INFINITY), _mm_cmple_ps(x, _mm_setzero_ps()));

    return p;
}


float dotlog_sse(const float* xs, const float* ys, const size_t n)
{
    __m128 xv, yv;
    union ans_t
    {
        __m128  v;
        float   f[4];
    } ans;
    ans.v = _mm_setzero_ps();

    size_t i;
    for (i = 0; i < n / 4; ++i) {
        xv = _mm_load_ps(xs + 4 * i);
        yv = _mm_load_ps(ys + 4 * i);
        ans.v = _mm_add_ps(ans.v, _mm_mul_ps(xv, log2_sse(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3];

    /* handle any overhang */
    i *= 4;
    switch (n % 4) {
        case 3: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 2: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 1: fans += xs[i] * fastlog2(ys[i]);
    }

    return fans;
}


float dotlogc_sse(const float* xs, const float* ys, const size_t n,
                  const float c)
{
    __m128 xv, yv;
    union ans_t
    {
        __m128  v;
        float   f[4];
    } ans;
    ans.v = _mm_setzero_ps();

    __m128 cv = _mm_set1_ps(c);

    size_t i;
    for (i = 0; i < n / 4; ++i) {
        xv = _mm_load_ps(xs + 4 * i);
        yv = _mm_mul_ps(cv, _mm_load_ps(ys + 4 * i));
        ans.v = _mm_add_ps(ans.v, _mm_mul_ps(xv, log2_sse(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3];

    /* handle any overhang */
    i *= 4;
    switch (n % 4) {
        case 3: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 2: fans += xs[i] * fastlog2(ys[i]); ++i;
        case 1: fans += xs[i] * fastlog2(ys[i]);
    }

    return fans;
}


void asxpy_sse(float* xs, const float* ys, const float c,
               const unsigned int* idx, const unsigned int off,
               const size_t n)
{
    __m128 yv;
    __m128 cv = _mm_set1_ps(c);
    union {
        __m128 v;
        float f[8];
    } x;

    union {
        __m128i v;
        unsigned int i[4];
    } iv;

    __m128i voff = _mm_set1_epi32((int) off);

    size_t i;
    for (i = 0; i < 4 * (n / 4); i += 4) {
        yv = _mm_mul_ps(cv, _mm_load_ps(ys + i));
        iv.v = _mm_load_si128(reinterpret_cast<const __m128i*>(idx + i));
        iv.v = _mm_sub_epi32(iv.v, voff);

        /* load from xs */
        x.f[0] = xs[iv.i[0]];
        x.f[1] = xs[iv.i[1]];
        x.f[2] = xs[iv.i[2]];
        x.f[3] = xs[iv.i[3]];

        x.v = _mm_add_ps(x.v, yv);

        /* store in xs */
        xs[iv.i[0]] = x.f[0];
        xs[iv.i[1]] = x.f[1];
        xs[iv.i[2]] = x.f[2];
        xs[iv.i[3]] = x.f[3];
    }

    /* handle overhang */
    i = 4 * (n / 4);
    switch (n % 4) {
        case 3: xs[idx[i] - off] += c * ys[i]; ++i;
        case 2: xs[idx[i] - off] += c * ys[i]; ++i;
        case 1: xs[idx[i] - off] += c * ys[i];
    }
}


float asxtydsz_sse(const float* xs, const float* ys, const float* zs,
                   const unsigned int* idx, const unsigned int off,
                   const size_t n)
{
    union {
        __m128 v;
        float f[4];
    } ans, x, z;
    ans.v = _mm_setzero_ps();

    union {
        __m128i a;
        unsigned int b[4];
    } i0;

    __m128i voff = _mm_set1_epi32((int) off);

    size_t i;
    for (i = 0; i < 4 * (n / 4); i += 4) {
        __m128 y  = _mm_load_ps(ys + i);

        i0.a = _mm_load_si128((__m128i*)(idx + i));
        i0.a = _mm_sub_epi32(i0.a, voff);

        x.f[0] = xs[i0.b[0]];
        x.f[1] = xs[i0.b[1]];
        x.f[2] = xs[i0.b[2]];
        x.f[3] = xs[i0.b[3]];

        ans.v = _mm_add_ps(ans.v, _mm_div_ps(_mm_mul_ps(x.v, y), z.v));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3];

    /* handle overhang */
    i = 4 * (n / 4);
    switch (n % 4) {
        case 3: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 2: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off]; ++i;
        case 1: fans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off];
    }

    return fans;
}


float dot_sse(const float* xs, const float* ys, const float* zs, size_t n)
{
    union ans_t
    {
        __m128  v;
        float   f[4];
    } ans;
    ans.v = _mm_setzero_ps();

    size_t i;
    __m128 x, y, z;
    for (i = 0; i < 4 * (n / 4); i += 4) {
        x = _mm_load_ps(xs + i);
        y = _mm_loadu_ps(ys + i);
        z = _mm_load_ps(zs + i);
        ans.v = _mm_add_ps(ans.v, _mm_mul_ps(_mm_mul_ps(x, y), z));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3];

    // handle overhang
    i = 4 * (n / 4);
    switch (n % 4) {
        case 3: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 2: fans += xs[i] * ys[i] * zs[i]; ++i;
        case 1: fans += xs[i] * ys[i] * zs[i];
    }

    return fans;
}


void* aalloc_vanilla(size_t n)
{
    void* xs = malloc(n);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree_vanilla(void* xs)
{
    free(xs);
}


void acopy_vanilla(void* dest, const void* src, size_t n)
{
    memcpy(dest, src, n);
}


float dotlog_vanilla(const float* xs, const float* ys, const size_t n)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(ys[i]);
    }
    return ans;
}


float dotlogc_vanilla(const float* xs, const float* ys, const size_t n, const float c)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(c * ys[i]);
    }
    return ans;

}


void asxpy_vanilla(float* xs, const float* ys, const float c,
                   const unsigned int* idx,
                   const unsigned int off,
                   const size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[idx[i] - off] += c * ys[i];
    }
}


float asxtydsz_vanilla(const float* xs, const float* ys, const float* zs,
                       const unsigned int* idx, const unsigned int off,
                       const size_t n)
{
    float ans = 0.0;
    for (size_t i = 0; i < n; ++i) {
        ans += xs[idx[i] - off] * ys[i] / zs[idx[i] - off];
    }
    return ans;
}


float dot_vanilla(const float* xs, const float* ys, const float* zs, const size_t n)
{
    float accum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        accum += xs[i] * ys[i] * zs[i];
    }
    return accum;
}


void linalg_init()
{
    if (cpu_has_avx()) {
        aalloc   = aalloc_avx;
        afree    = afree_avx;
        acopy    = acopy_avx;
        dotlog   = dotlog_avx;
        dotlogc  = dotlogc_avx;
        asxpy    = asxpy_avx;
        asxtydsz = asxtydsz_avx;
        dot      = dot_avx;
        LINALG_INSTR_SET = "AVX";
    }
    else if (cpu_has_sse2()) {
        aalloc   = aalloc_sse;
        afree    = afree_sse;
        acopy    = acopy_sse;
        dotlog   = dotlog_sse;
        dotlogc  = dotlogc_sse;
        asxpy    = asxpy_sse;
        asxtydsz = asxtydsz_sse;
        dot      = dot_sse;
        if (cpu_has_sse4()) {
            log2_sse = log2_sse4;
            LINALG_INSTR_SET = "SSE4";
        }
        else {
            log2_sse = log2_sse2;
            LINALG_INSTR_SET = "SSE2";
        }
    }
    else {
        aalloc   = aalloc_vanilla;
        afree    = afree_vanilla;
        acopy    = acopy_vanilla;
        dotlog   = dotlog_vanilla;
        dotlogc  = dotlogc_vanilla;
        asxpy    = asxpy_vanilla;
        asxtydsz = asxtydsz_vanilla;
        dot      = dot_vanilla;
        LINALG_INSTR_SET = "Vanilla";
    }
}



