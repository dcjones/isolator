
#include <cmath>

#include "constants.hpp"
#include "fastmath.hpp"
#include "fastmath-vanilla.hpp"
#include "logger.hpp"


static const float prob_epsilon = constants::frag_prob_epsilon;


void* aalloc_vanilla(size_t n)
{
    void* xs = malloc(n);
    if (xs == NULL) {
        Logger::abort("Can't allocate an array of size %ul.",
                      (unsigned long) n);
    }
    return xs;
}


void afree_vanilla(void* xs)
{
    free(xs);
}


void acopy_vanilla(void* dest, const void* src, size_t n)
{
    memcpy(dest, src, n);
}


float sumlog_vanilla(const float* xs, const size_t n)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += fastlog2(xs[i]);
    }
    return ans;
}


float dotlog_vanilla(const float* xs, const float* ys, const size_t n)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(ys[i]);
    }
    return ans;
}


float dotlogc_vanilla(const float* xs, const float* ys, const size_t n, const float c)
{
    float ans = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        ans += xs[i] * fastlog2(c * ys[i]);
    }
    return ans;

}


void asxpy_vanilla(float* xs, const float* ys, const float c,
                   const unsigned int* idx,
                   const unsigned int off,
                   const size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[idx[i] - off] =
            std::max<float>(xs[idx[i] - off] + c * ys[i], prob_epsilon);
    }
}


void axpy_vanilla(float* xs, const float* ys, const float c, const size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        xs[i] = std::max<float>(xs[i] + c * ys[i], prob_epsilon);
    }
}


float asxtydsz_vanilla(const float* ys, const float* zs,
                       const unsigned int* idx, const unsigned int off,
                       const size_t n)
{
    float ans = 0.0;
    for (size_t i = 0; i < n; ++i) {
        ans += ys[i] / zs[idx[i] - off];
    }
    return ans;
}


float sumdiv_vanilla(const float* xs, const float* ys, const size_t n)
{
    float ans = 0.0;
    for (size_t i = 0; i < n; ++i) {
        ans += xs[i] / ys[i];
    }
    return ans;
}


float dot_vanilla(const float* xs, const float* ys, const float* zs, const size_t n)
{
    float accum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        accum += xs[i] * ys[i] * zs[i];
    }
    return accum;
}


