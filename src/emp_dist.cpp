
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>

#include "common.hpp"
#include "emp_dist.hpp"


EmpDist::EmpDist(EmpDist& other)
    : n(other.n)
    , m(other.m)
    , med(other.med)
    , pdf_memo(other.pdf_memo)
    , cdf_memo(other.cdf_memo)
{
    vals = new unsigned int [n];
    lens = new unsigned int [n];

    std::copy(other.vals, other.vals + n, vals);
    std::copy(other.lens, other.lens + n, lens);
}


EmpDist::EmpDist(const unsigned int* vals, const unsigned int* lens, size_t n)
    : med(NAN)
{
    /* count number of non-zero observations. */
    size_t nnz = 0;
    for (size_t i = 0; i < n; ++i) {
        if (lens[i] > 0) nnz += 1;
    }

    this->vals = new unsigned int [nnz];
    this->lens = new unsigned int [nnz];

    m = 0;
    for (size_t i = 0, j = 0; i < n; ++i) {
        if (lens[i] > 0) {
            this->vals[j] = vals[i];
            this->lens[j] = lens[i];
            m += lens[i];
            ++j;
        }
    }

    this->n = nnz;
    this->med = compute_median();
}


float EmpDist::compute_median() const
{
    if (n == 0) return NAN;

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

    return vals[i];
}


EmpDist::~EmpDist()
{
    delete [] vals;
    delete [] lens;
}


float EmpDist::median() const
{
    if (n == 0) return NAN;
    else return med;
}


float EmpDist::eval(unsigned int x) const
{
    if (x > vals[n - 1]) return 0.0;

    boost::unordered_map<unsigned int, float>::iterator idx;
    idx = pdf_memo.find(x);
    if (idx != pdf_memo.end()) {
        return idx->second;
    }

    /* binary search for the nearest item in vals */
    unsigned int i = std::lower_bound(vals, vals + n, x) - vals;
    float ans = vals[i] == x ? lens[i] / (float) m : 0.0;

    pdf_memo.insert(std::pair<unsigned int, double>(x, ans));
    return ans;
}


float EmpDist::pdf(unsigned int x) const
{
    boost::lock_guard<boost::mutex> lock(mut);
    return eval(x);
}


float EmpDist::cdf(unsigned int x) const
{
    boost::lock_guard<boost::mutex> lock(mut);

    boost::unordered_map<unsigned int, float>::iterator idx;
    idx = cdf_memo.find(x);
    if (idx != cdf_memo.end()) return idx->second;

    float ans = 0.0;
    size_t i;
    for (i = 0; i < n && vals[i] <= x; ++i) {
        ans += (float) lens[i];
    }
    ans /= (float) m;

    cdf_memo.insert(std::make_pair(x, ans));

    return ans;
}


