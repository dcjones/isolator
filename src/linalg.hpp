
#ifndef ISOLATOR_LINALG_HPP
#define ISOLATOR_LINALG_HPP

#include <cstdlib>

/* A couple very specific linear algebra functions implemented using
 * different instruction sets (vanilla x86, SSE4, and AVX). */


/* Allocate a vector, with proper alignment if required. */
float* vector_alloc(size_t n);
void vector_free(float*);


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
 *     xs[idx[i]] += c * ys[i]
 *
 *  For all 0 <= i < n.
 */
void asxpy(float* xs, const float* ys, const float c,
            const unsigned int* idx, const size_t n);


#endif

