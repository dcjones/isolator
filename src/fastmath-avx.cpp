
#include "config.h"

#include <cmath>
#include <immintrin.h>

#include "constants.hpp"
#include "fastmath.hpp"
#include "logger.hpp"


#ifdef _MSC_VER
  #define ALIGN16_START __declspec(align(16))
  #define ALIGN16_END
  #define ALIGN32_START __declspec(align(32))
  #define ALIGN32_END
#else
  #define ALIGN16_START
  #define ALIGN16_END __attribute__((aligned(16)))
  #define ALIGN32_START
  #define ALIGN32_END __attribute__((aligned(32)))
#endif

#define PS16_CONST(name, c) \
    static const ALIGN16_START float ps16_##name[8] ALIGN16_END = {c, c, c, c}

#define PI16_CONST(name, c) \
    static const ALIGN16_START int pi16_##name[8] ALIGN16_END = {c, c, c, c}

#define PS32_CONST(name, c) \
    static const ALIGN32_START float ps32_##name[8] ALIGN32_END = {c, c, c, c, c, c, c, c}

#define PI32_CONST(name, c) \
    static const ALIGN32_START int pi32_##name[8] ALIGN32_END = {c, c, c, c, c, c, c, c}


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
    for (i = 0; i + 8 <= n; i += 8) {
        x = _mm256_load_ps(src + i);
        _mm256_store_ps(dest + i, x);
    }

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
    PS32_CONST(1, 1.0f);
    PI32_CONST(inv_mant_mask, ~0x7f800000);

    PS32_CONST(log2_c0, 3.1157899f);
    PS32_CONST(log2_c1, -3.3241990f);
    PS32_CONST(log2_c2, 2.5988452f);
    PS32_CONST(log2_c3, -1.2315303f);
    PS32_CONST(log2_c4, 3.1821337e-1f);
    PS32_CONST(log2_c5, -3.4436006e-2f);
    PS32_CONST(neginf, (float) -INFINITY);

    PS32_CONST(exp_c, 1.1920928955078125e-7f);
    PI32_CONST(mant_mask, 0x7f800000);
    PS32_CONST(7f, 127.0f);

    /* extract the exponent */
    __m256 e = x;
    e = _mm256_and_ps(e, *(__m256*) pi32_mant_mask);
    e = _mm256_cvtepi32_ps( *(__m256i*) &e);
    e = _mm256_mul_ps(e, *(__m256*) ps32_exp_c);
    e = _mm256_sub_ps(e, *(__m256*) ps32_7f);

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
    for (i = 0; i + 8 <= n; i += 8) {
        xv = _mm256_load_ps(xs + i);
        yv = _mm256_load_ps(ys + i);
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(xv, log2_avx(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3] +
                 ans.f[4] + ans.f[5] + ans.f[6] + ans.f[7];

    /* handle any overhang */
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
    static const float prob_epsilon = constants::frag_prob_epsilon;
    PS32_CONST(prob_epsilon, prob_epsilon);

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
    for (i = 0; i + 8 <= n; i += 8) {
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
        x.v = _mm256_max_ps(x.v, *reinterpret_cast<const __m256*>(ps32_prob_epsilon));

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
    switch (n % 8) {
        case 7:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 6:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 5:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 4:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 3:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 2:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
            ++i;
        case 1:
            xs[idx[i] - off] = std::max<float>(prob_epsilon, xs[idx[i] - off] + c * ys[i]);
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
    for (i = 0; i + 8 <= n; i += 8) {
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


