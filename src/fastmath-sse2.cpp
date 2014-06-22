

#include "config.h"

#include <cmath>
#include <immintrin.h>

#include "constants.hpp"
#include "fastmath.hpp"
#include "logger.hpp"
#include "cpuid.hpp"


#ifdef _MSC_VER
  #define ALIGN16_START __declspec(align(16))
  #define ALIGN16_END
#else
  #define ALIGN16_START
  #define ALIGN16_END __attribute__((aligned(16)))
#endif

#define PS16_CONST(name, c) \
    static const ALIGN16_START float ps16_##name[8] ALIGN16_END = {c, c, c, c}

#define PI16_CONST(name, c) \
    static const ALIGN16_START int pi16_##name[8] ALIGN16_END = {c, c, c, c}

// set in fastmath_sse_init
static __m128 (*log2_sse)(__m128 x) = NULL;

void* aalloc_sse(size_t n)
{
    void* xs = _mm_malloc(n, 32);
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


/* Comptue log2 over an sse single vector. */
__m128 log2_sse2(__m128 x)
{
    PS16_CONST(1, 1.0f);
    PI16_CONST(inv_mant_mask, ~0x7f800000);
    PS16_CONST(log2_c0, 3.1157899f);
    PS16_CONST(log2_c1, -3.3241990f);
    PS16_CONST(log2_c2, 2.5988452f);
    PS16_CONST(log2_c3, -1.2315303f);
    PS16_CONST(log2_c4, 3.1821337e-1f);
    PS16_CONST(log2_c5, -3.4436006e-2f);
    PS16_CONST(neginf, (float) -INFINITY);

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
    // TODO: test this
    __m128 mask = _mm_cmpgt_ps(x, _mm_setzero_ps());
    p = _mm_or_ps(_mm_and_ps(mask, p), _mm_andnot_ps(mask, *(__m128*) ps16_neginf));

    return p;
}


// defined in fastmath-sse4.cpp
__m128 log2_sse4(__m128 x);


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
    static const float prob_epsilon = constants::frag_prob_epsilon;
    PS16_CONST(prob_epsilon, prob_epsilon);

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
        x.v = _mm_max_ps(x.v, *reinterpret_cast<const __m128*>(ps16_prob_epsilon));

        /* store in xs */
        xs[iv.i[0]] = x.f[0];
        xs[iv.i[1]] = x.f[1];
        xs[iv.i[2]] = x.f[2];
        xs[iv.i[3]] = x.f[3];
    }

    /* handle overhang */
    i = 4 * (n / 4);
    switch (n % 4) {
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

        z.f[0] = zs[i0.b[0]];
        z.f[1] = zs[i0.b[1]];
        z.f[2] = zs[i0.b[2]];
        z.f[3] = zs[i0.b[3]];

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


void fastmath_sse_init()
{
    if (cpu_has_sse4()) {
        log2_sse = log2_sse4;
    }
    else {
        log2_sse = log2_sse2;
    }
}



