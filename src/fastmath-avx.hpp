
#ifndef ISOLATOR_FASTMATH_AVX_HPP
#define ISOLATOR_FASTMATH_AVX_HPP

#include <cstdlib>

void* aalloc_avx(size_t n);
void afree_avx(void* xs);
void acopy_avx(void* dest, const void* src, size_t n);
float dotlog_avx(const float* xs, const float* ys, const size_t n);
float dotlogc_avx(const float* xs, const float* ys, const size_t n, const float c);
void asxpy_avx(float* xs, const float* ys, const float c,
               const unsigned int* idx,
               const unsigned int off,
               const size_t n);
float asxtydsz_avx(const float* xs, const float* ys, const float* zs,
                   const unsigned int* idx, const unsigned int off,
                   const size_t n);
float dot_avx(const float* xs, const float* ys, const float* zs, const size_t n);

#endif

