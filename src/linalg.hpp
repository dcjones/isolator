
#ifndef ISOLATOR_LINALG_HPP
#define ISOLATOR_LINALG_HPP

#include <cstdlib>

/* A couple very specific linear algebra functions implemented using
 * different instruction sets (vanilla x86, SSE4, and AVX). */


/* Allocate an array with 16-bytes alignment (assuming the cpu supports SEE or
 * AVX) */
void* aalloc(size_t n);
void afree(void*);

/* Fast copying of aligned arrays.
 * Probably most memcpy implementations will do this sort of optimization
 * automatically, but being able to assume alignment might buy us a little.
 * */
void acopy(void* dest, const void* src, size_t n);

/* Copy an array of n aligned doubles to an array of n aligned floats. */
void pdpscpy(float* dest, const double* src, size_t n);

/* Copy an array of n aligned floats to an array of n aligned doubles. */
void pspscpy(float* dest, const double* src, size_t n);


/* Dot product of xs and log(ys).
 *
 * That is,
 *  xs[0] * log(ys[0]) + ... + xs[n - 1] * log(ys[n - 1])
 *
 * */
float dotlog(const float* xs, const float* ys, const size_t n);


/* Weighted sparse vector addition.
 *
 * If xs is a dense vector and ys a sparse vector with n non-zero entries at
 * indexes idx[0], ..., idx[n - 1], compute:
 *
 *     xs[idx[i] - off] += c * ys[i]
 *
 *  For all 0 <= i < n.
 */
void asxpy(float* xs, const float* ys, const float c,
           const unsigned int* idx, const unsigned int off,
           const size_t n);


#endif

