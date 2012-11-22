
#include "config.h"
#include "linalg.hpp"
#include "logger.hpp"

#include <math.h>

#define POLY0(x, c0) (c0)
#define POLY1(x, c0, c1) (x * POLY0(x, c1) + c0)
#define POLY2(x, c0, c1, c2) (x * POLY1(x, c1, c2) + c0)
#define POLY3(x, c0, c1, c2, c3) (x * POLY2(x, c1, c2, c3) + c0)
#define POLY4(x, c0, c1, c2, c3, c4) (x * POLY3(x, c1, c2, c3, c4) + c0)
#define POLY5(x, c0, c1, c2, c3, c4, c5) (x * POLY4(x, c1, c2, c3, c4, c5) + c0)

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

    float p = POLY5(m, 3.1157899f, -3.3241990f, 2.5988452f, -1.2315303f,  3.1821337e-1f, -3.4436006e-2f);

    p *= m - 1.0;

    return p + e;
}



/* AVX versions */
#if defined(HAVE_AVX) && defined(HAVE_IMMINTRIN_H) && defined(__AVX__)

#include <immintrin.h>


void* aalloc(size_t n)
{
    void* xs = _mm_malloc(n, 32);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree(void* xs)
{
    _mm_free(xs);
}


void acopy(void* dest_, const void* src_, size_t n)
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


typedef union {
    __m256i a;
    __m128i b[2];
} m256i_m128i_t;

typedef union {
    __m256 a;
    __m128 b[2];
} m256_m128_t;

/* Macros for evaluating ploynomials. */
#define AVX_POLY0(x, c0) _mm256_set1_ps(c0)
#define AVX_POLY1(x, c0, c1) \
    _mm256_add_ps(_mm256_mul_ps(AVX_POLY0(x, c1), x), _mm256_set1_ps(c0))
#define AVX_POLY2(x, c0, c1, c2) \
    _mm256_add_ps(_mm256_mul_ps(AVX_POLY1(x, c1, c2), x), _mm256_set1_ps(c0))
#define AVX_POLY3(x, c0, c1, c2, c3) \
    _mm256_add_ps(_mm256_mul_ps(AVX_POLY2(x, c1, c2, c3), x), _mm256_set1_ps(c0))
#define AVX_POLY4(x, c0, c1, c2, c3, c4) \
    _mm256_add_ps(_mm256_mul_ps(AVX_POLY3(x, c1, c2, c3, c4), x), _mm256_set1_ps(c0))
#define AVX_POLY5(x, c0, c1, c2, c3, c4, c5) \
    _mm256_add_ps(_mm256_mul_ps(AVX_POLY4(x, c1, c2, c3, c4, c5), x), _mm256_set1_ps(c0))


/* Comptue log2 over an avx single vector. */
static __m256 avx_log2(__m256 x)
{
    /* TODO: this can be written using avx2 instructions, without resorting to
     * SSE for integer operations, but I won't bother until I actually have an
     * cpu I can test that on.
     * */

    /* extract the exponent */
    const __m128 c127 = _mm_set1_epi32(127);
    m256i_m128i_t i;
    m256_m128_t e;
    i.a = _mm256_castps_si256(x);
    e.b[0] = _mm_sub_epi32(_mm_srli_epi32(i.b[0], 23), c127);
    e.b[1] = _mm_sub_epi32(_mm_srli_epi32(i.b[1], 23), c127);
    e.a = _mm256_cvtepi32_ps(e.a);

    /* extract the mantissa */
    m256_m128_t m;
    const __m256 c1 = _mm256_set1_ps(1.0f);
    const __m128i mant_mask = _mm_set1_epi32(0x007FFFFF);
    m.b[0] = _mm_castsi128_ps(_mm_and_si128(i.b[0], mant_mask));
    m.b[1] = _mm_castsi128_ps(_mm_and_si128(i.b[1], mant_mask));
    m.a = _mm256_or_ps(m.a, c1);

    /* polynomial approximation on the mantissa */
    __m256 p = AVX_POLY5(m.a, 3.1157899f,
                             -3.3241990f,
                              2.5988452f,
                             -1.2315303f,
                              3.1821337e-1f,
                             -3.4436006e-2f);

    p = _mm256_add_ps(_mm256_mul_ps(p, _mm256_sub_ps(m.a, c1)), e.a);

    /* Handle the log(x) = -INFINITY, for x <= 0 cases. */
    p = _mm256_blendv_ps(p, _mm256_set1_ps(-INFINITY),
                         _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LE_OQ));

    return p;
}


float dotlog(const float* xs, const float* ys, const size_t n)
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
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(xv, avx_log2(yv)));
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


float dotlogc(const float* xs, const float* ys, const size_t n, const float c)
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
        ans.v = _mm256_add_ps(ans.v, _mm256_mul_ps(xv, avx_log2(yv)));
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


void asxpy(float* xs, const float* ys, const float c,
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
        __m256i v;
        __m128i w[2];
        unsigned int i[8];
    } iv;

    __m128i voff = _mm_set1_epi32((int) off);

    size_t i;
    for (i = 0; i < 8 * (n / 8); i += 8) {
        yv = _mm256_mul_ps(cv, _mm256_load_ps(ys + i));
        iv.v = _mm256_load_si256(reinterpret_cast<const __m256i*>(idx + i));
        iv.w[0] = _mm_sub_epi32(iv.w[0], voff);
        iv.w[1] = _mm_sub_epi32(iv.w[1], voff);

        /* load from xs */
        x.f[0] = xs[iv.i[0]];
        x.f[1] = xs[iv.i[1]];
        x.f[2] = xs[iv.i[2]];
        x.f[3] = xs[iv.i[3]];
        x.f[4] = xs[iv.i[4]];
        x.f[5] = xs[iv.i[5]];
        x.f[6] = xs[iv.i[6]];
        x.f[7] = xs[iv.i[7]];

        x.v = _mm256_add_ps(x.v, yv);

        /* store in xs */
        xs[iv.i[0]] = x.f[0];
        xs[iv.i[1]] = x.f[1];
        xs[iv.i[2]] = x.f[2];
        xs[iv.i[3]] = x.f[3];
        xs[iv.i[4]] = x.f[4];
        xs[iv.i[5]] = x.f[5];
        xs[iv.i[6]] = x.f[6];
        xs[iv.i[7]] = x.f[7];
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


/* SSE2 versions */
#elif defined(HAVE_IMMINTRIN_H) && defined(HAVE_SSE2) && defined(__SSE2__)

#include <immintrin.h>

void* aalloc(size_t n)
{
    void* xs = _mm_malloc(n, 16);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree(void* xs)
{
    _mm_free(xs);
}


void acopy(void* dest_, const void* src_, size_t n)
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


/* Macros for evaluating ploynomials. */
#define SSE_POLY0(x, c0) _mm_set1_ps(c0)
#define SSE_POLY1(x, c0, c1) \
    _mm_add_ps(_mm_mul_ps(SSE_POLY0(x, c1), x), _mm_set1_ps(c0))
#define SSE_POLY2(x, c0, c1, c2) \
    _mm_add_ps(_mm_mul_ps(SSE_POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define SSE_POLY3(x, c0, c1, c2, c3) \
    _mm_add_ps(_mm_mul_ps(SSE_POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define SSE_POLY4(x, c0, c1, c2, c3, c4) \
    _mm_add_ps(_mm_mul_ps(SSE_POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define SSE_POLY5(x, c0, c1, c2, c3, c4, c5) \
    _mm_add_ps(_mm_mul_ps(SSE_POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

/* Comptue log2 over an sse single vector. */
static __m128 sse_log2(__m128 x)
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
    __m128 p = SSE_POLY5(m,
                         3.1157899f,
                        -3.3241990f,
                         2.5988452f,
                        -1.2315303f,
                         3.1821337e-1f,
                        -3.4436006e-2f);

    p = _mm_add_ps(_mm_mul_ps(p, _mm_sub_ps(m, one)), e);

    /* handle the case with x= < 0.0 */
#if (defined(HAVE_SSE41) && defined(__SSE41__)) || (defined(HAVE_SSE42) && defined(__SSE42__))
    p = _mm_blendv_ps(p, _mm_set1_ps(-INFINITY), _mm_cmple_ps(x, _mm_setzero_ps()));
#else
    // TODO: test this
    __m128 mask = _mm_cmpgt_ps(x, _mm_setzero_ps());
    p = _mm_or_ps(_mm_and_ps(mask, p), _mm_andnot_ps(mask, _mm_set1_ps(-INFINITY)));
#endif

    return p;
}


float dotlog(const float* xs, const float* ys, const size_t n)
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
        ans.v = _mm_add_ps(ans.v, _mm_mul_ps(xv, sse_log2(yv)));
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


float dotlog(const float* xs, const float* ys, const size_t n, const float c)
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
        ans.v = _mm_add_ps(ans.v, _mm_mul_ps(xv, sse_log2(yv)));
    }

    float fans = ans.f[0] + ans.f[1] + ans.f[2] + ans.f[3];

    /* handle any overhang */
    i *= 4;
    switch (n % 4) {
        case 3: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 2: fans += xs[i] * fastlog2(c * ys[i]); ++i;
        case 1: fans += xs[i] * fastlog2(c * ys[i]);
    }

    return fans;
}


void asxpy(float* xs, const float* ys, const float c,
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


/* Vanilla versions */
#else


void* aalloc(size_t n)
{
    void* xs = malloc(n);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree(void* xs)
{
    free(xs);
}


void acopy(void* dest, const void* src, size_t n)
{
    memcpy(dest, src, n);
}


float dotlog(const float* xs, const float* ys, const size_t n)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(ys[i]);
    }
    return ans;
}


float dotlogc(const float* xs, const float* ys, const size_t n, const float c)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(c * ys[i]);
    }
    return ans;

}


void asxpy(float* xs, const float* ys, const float c,
            const unsigned int* idx,
            const unsigned int off,
            const size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[idx[i] - off] += c * ys[i];
    }
}


#endif
