
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


/* Dot product of xs and log(ys).
 *
 * That is,
 *  xs[0] * log(ys[0]) + ... + xs[n - 1] * log(ys[n - 1])
 *
 * */
float dotlog(const float* xs, const float* ys, const size_t n);


/* Dot product of xs and log(c * ys).
 * That is,
 *  xs[0] * log(c * ys[0]) + ... + xs[n - 1] * log(c * ys[n - 1])
 */
float dotlogc(const float* xs, const float* ys, const size_t n, const float c);


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

/* ssxtydsz: sum of sparse x times y divided by sparse z.
 *
 * ...I know, I know, but ou try coming up with a better name!
 *
 * This computes the sum of
 *   xs[idx[i] - off] * ys[i] / zs[idx[i] - off]
 *
 * for 0 <= i < n.
 */
float asxtydsz(const float* xs, const float* ys, const float* zs,
               const unsigned int* idx, const unsigned int off,
               const size_t n);


/* Fast log2 approximation. */
float fastlog2(float x);

#endif

