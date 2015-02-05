
#include <cmath>

#include "fastmath.hpp"
#include "fastmath-vanilla.hpp"
#include "fastmath-sse.hpp"
#include "fastmath-avx.hpp"
#include "cpuid.hpp"

/* The implementation of each of these functions is determined at run time */
void* (*aalloc)(size_t n) = NULL;
void (*afree)(void*) = NULL;
void (*acopy)(void* dest, const void* src, size_t n) = NULL;
float (*sumlog)(const float* xs, const size_t n) = NULL;
float (*dotlog)(const float* xs, const float* ys, const size_t n) = NULL;
float (*dotlogc)(const float* xs, const float* ys, const size_t n, const float c) = NULL;
void (*asxpy)(float* xs, const float* ys, const float c,
              const unsigned int* idx, const unsigned int off,
              const size_t n) = NULL;
void (*axpy)(float* xs, const float* ys, const float c, const size_t n) = NULL;
float (*asxtydsz)(const float* ys, const float* zs,
                  const unsigned int* idx, const unsigned int off,
                  const size_t n) = NULL;
float (*sumdiv)(const float* xs, const float* ys, const size_t n);
float (*dot)(const float* xs, const float* ys, const float* zs, size_t n);
const char* FASTMATH_INSTR_SET = "";


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


void fastmath_init()
{
    if (cpu_has_avx()) {
        aalloc   = aalloc_avx;
        afree    = afree_avx;
        acopy    = acopy_avx;
        sumlog   = sumlog_avx;
        dotlog   = dotlog_avx;
        dotlogc  = dotlogc_avx;
        asxpy    = asxpy_avx;
        axpy     = axpy_avx;
        asxtydsz = asxtydsz_avx;
        sumdiv   = sumdiv_avx;
        dot      = dot_avx;
        FASTMATH_INSTR_SET = "AVX";
    }
    else if (cpu_has_sse2()) {
        aalloc   = aalloc_sse;
        afree    = afree_sse;
        acopy    = acopy_sse;
        sumlog   = sumlog_sse;
        dotlog   = dotlog_sse;
        dotlogc  = dotlogc_sse;
        asxpy    = asxpy_sse;
        axpy     = axpy_vanilla; // TODO
        asxtydsz = asxtydsz_sse;
        sumdiv   = sumdiv_vanilla; // TODO
        dot      = dot_sse;
        fastmath_sse_init();
        if (cpu_has_sse4()) {
            FASTMATH_INSTR_SET = "SSE4";
        }
        else {
            FASTMATH_INSTR_SET = "SSE2";
        }
    }
    else {
        aalloc   = aalloc_vanilla;
        afree    = afree_vanilla;
        acopy    = acopy_vanilla;
        sumlog   = sumlog_vanilla;
        dotlog   = dotlog_vanilla;
        dotlogc  = dotlogc_vanilla;
        asxpy    = asxpy_vanilla;
        axpy     = axpy_vanilla;
        asxtydsz = asxtydsz_vanilla;
        sumdiv   = sumdiv_vanilla;
        dot      = dot_vanilla;
        FASTMATH_INSTR_SET = "Vanilla";
    }
}


