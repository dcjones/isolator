
#ifndef ISOLATOR_FASTMATH_VANILLA_HPP
#define ISOLATOR_FASTMATH_VANILLA_HPP

#include <cstdlib>

void* aalloc_vanilla(size_t n);
void afree_vanilla(void* xs);
void acopy_vanilla(void* dest, const void* src, size_t n);
float dotlog_vanilla(const float* xs, const float* ys, const size_t n);
float dotlogc_vanilla(const float* xs, const float* ys, const size_t n, const float c);
void asxpy_vanilla(float* xs, const float* ys, const float c,
                   const unsigned int* idx,
                   const unsigned int off,
                   const size_t n);
float asxtydsz_vanilla(const float* xs, const float* ys, const float* zs,
                       const unsigned int* idx, const unsigned int off,
                       const size_t n);
float dot_vanilla(const float* xs, const float* ys, const float* zs, const size_t n);

#endif


