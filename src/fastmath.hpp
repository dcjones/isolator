
#ifndef ISOLATOR_FASTMATH_HPP
#define ISOLATOR_FASTMATH_HPP

#include <cstdlib>

/* A few very specific linear algebra functions implemented using
 * different instruction sets (vanilla x86, SSE4, and AVX). */


/* This function must be called before any of the functions below are. It is
 * responsible for setting the otherwise NULL function pointers. */
void fastmath_init();


/* String describing the instruction set used. Set after calling fastmath_init. */
extern const char* FASTMATH_INSTR_SET;


/* Allocate an array with 16-bytes alignment (assuming the cpu supports SEE or
 * AVX) */
extern void* (*aalloc)(size_t n);
extern void (*afree)(void*);

/* Fast copying of aligned arrays.
 * Probably most memcpy implementations will do this sort of optimization
 * automatically, but being able to assume alignment might buy us a little.
 * */
extern void (*acopy)(void* dest, const void* src, size_t n);


/* Dot product of xs and log(ys).
 *
 * That is,
 *  xs[0] * log(ys[0]) + ... + xs[n - 1] * log(ys[n - 1])
 *
 * */
extern float (*dotlog)(const float* xs, const float* ys, const size_t n);


/* Sum of logorithms of xs.
 *
 * That is,
 *  log(xs[0]) + ... + log(xs[n - 1])
 */
extern float (*sumlog)(const float* xs, const size_t n);


/* Dot product of xs and log(c * ys).
 * That is,
 *  xs[0] * log(c * ys[0]) + ... + xs[n - 1] * log(c * ys[n - 1])
 */
extern float (*dotlogc)(const float* xs, const float* ys, const size_t n, const float c);


/* Weighted sparse vector addition.
 *
 * If xs is a dense vector and ys a sparse vector with n non-zero entries at
 * indexes idx[0], ..., idx[n - 1], compute:
 *
 *     xs[idx[i] - off] += c * ys[i]
 *
 *  For all 0 <= i < n.
 */
extern void (*asxpy)(float* xs, const float* ys, const float c,
                     const unsigned int* idx, const unsigned int off,
                     const size_t n);

/*
 * Compute.
 *    xs[i] += c * ys[i]
 *
 */
extern void (*axpy)(float* xs, const float* ys, const float c, const size_t n);


/* ssxtydsz: sum of sparse x times y divided by sparse z.
 *
 * ...I know, I know, but you try coming up with a better name!
 *
 * This computes the sum of
 *   ys[i] / zs[idx[i] - off]
 *
 * for 0 <= i < n.
 */
extern float (*asxtydsz)(const float* ys, const float* zs,
                         const unsigned int* idx, const unsigned int off,
                         const size_t n);

/*
 * Compute the sum of:
 *   xs[i] / ys[i]
 */
extern float (*sumdiv)(const float* xs, const float* ys, const size_t n);


/* Triple product, where xs and zs are aligned and ys may not be.
 */
 extern float (*dot)(const float* xs, const float* ys, const float* zs, const size_t n);


/* Fast log2 approximation. */
float fastlog2(float x);
float fastlog(float x);

#endif

