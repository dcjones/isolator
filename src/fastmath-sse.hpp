
#ifndef ISOLATOR_FASTMATH_SSE_HPP
#define ISOLATOR_FASTMATH_SSE_HPP

#include <cstdlib>

void* aalloc_sse(size_t n);
void afree_sse(void* xs);
void acopy_sse(void* dest, const void* src, size_t n);
float dotlog_sse(const float* xs, const float* ys, const size_t n);
float dotlogc_sse(const float* xs, const float* ys, const size_t n, const float c);
void asxpy_sse(float* xs, const float* ys, const float c,
               const unsigned int* idx,
               const unsigned int off,
               const size_t n);
float asxtydsz_sse(const float* xs, const float* ys, const float* zs,
                   const unsigned int* idx, const unsigned int off,
                   const size_t n);
float dot_sse(const float* xs, const float* ys, const float* zs, const size_t n);
void fastmath_sse_init();

#endif


