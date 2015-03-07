
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>

#include "common.hpp"
#include "emp_dist.hpp"


EmpDist::EmpDist(EmpDist& other)
    : pdfvals(other.pdfvals)
    , cdfvals(other.cdfvals)
    , med(other.med)
{
}


EmpDist::EmpDist(const unsigned int* vals, const unsigned int* lens, size_t n)
{
    double valsum = 0;
    for (size_t i = 0; i < n; ++i) {
        valsum += lens[i];
    }

    double cumpr = 0.0;
    unsigned int lastval = 0, maxval = 1;
    for (; lastval < n; ++lastval) {
        cumpr += lens[lastval] / valsum;
        maxval = vals[lastval];
        if (cumpr > 1.0 - 1e-6) {
            break;
        }
    }

    pdfvals.resize(maxval);
    valsum = 0.0;
    for (unsigned int i = 0; i < lastval; ++i) {
        valsum += lens[i];
    }

    for (unsigned int val = 0, i = 0; val < maxval; ) {
        if (val == vals[i]) {
            pdfvals[val] = lens[i] / valsum;
            ++val;
            ++i;
        }
        else if (val < vals[i]) {
            pdfvals[val] = 0.0;
            ++val;
        }
    }

    cdfvals.resize(maxval);
    cdfvals[0] = pdfvals[0];
    for (unsigned int val = 1; val < maxval; ++val) {
        cdfvals[val] = cdfvals[val - 1] + pdfvals[val];
    }

    // compute median
    size_t i = 0, j = n - 1;
    unsigned int u = lens[0], v = lens[n - 1];
    while (i < j) {
        if (u <= v) {
            v -= u;
            u = lens[++i];
        }
        else {
            u -= v;
            v = lens[--j];
        }
    }
    med = vals[i];
}


float EmpDist::median() const
{
    if (pdfvals.size() == 0) return NAN;
    else return med;
}


float EmpDist::pdf(unsigned int x) const
{
    return x < pdfvals.size() ? pdfvals[x] : 0.0;
}


float EmpDist::cdf(unsigned int x) const
{
    return x < cdfvals.size() ? cdfvals[x] : 1.0;
}


