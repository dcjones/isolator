
#include "config.h"

#include <cmath>
#include <immintrin.h>


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


__m128 log2_sse4(__m128 x)
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
    p = _mm_blendv_ps(p, *(__m128*) ps16_neginf, _mm_cmple_ps(x, _mm_setzero_ps()));

    return p;
}


