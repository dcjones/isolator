
#include "common.hpp"
#include "emp_dist.hpp"


EmpDist::EmpDist(EmpDist& other)
    : n(other.n)
    , m(other.m)
    , w(other.w)
    , med(other.med)
    , pdf_memo(other.pdf_memo)
    , cdf_memo(other.cdf_memo)
{
    vals = new unsigned int [n];
    lens = new unsigned int [n];

    std::copy(other.vals, other.vals + n, vals);
    std::copy(other.lens, other.lens + n, lens);
}


EmpDist::EmpDist(const unsigned int* vals, const unsigned int* lens,
                 size_t n, float w)
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
    this->w = w;
}


EmpDist::~EmpDist()
{
    delete vals;
    delete lens;
}


float EmpDist::median() const
{
    if (n == 0) return NAN;
    if (finite(med)) return med;

    /* Dumb O(n) median computation. */
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

    return med = (float) vals[i];
}


/* A geometric weight function, as in,
 *
 * Wang, M. C., & Van Ryzin, J. (1981). A class of smooth estimators for
 * discrete distributions. Biometrika, 68(1), 301. Biometrika Trust.
 */
static float wf_geom(unsigned int u, unsigned int v, float w)
{
    if (u == v) {
        return 1.0 - w;
    }
    else {
        float x = abs((float) u - (float) v);
        return 0.5 * (1.0 - w) * powf(w, x);
    }
}


float EmpDist::eval(unsigned int x) const
{
    /* An approximation made for efficiency. Assign a probability of 0 to a
     * value that is greater than any of our observations. */
    if (x > vals[n - 1]) return 0.0;

    boost::unordered_map<unsigned int, float>::iterator idx;
    idx = pdf_memo.find(x);
    if (idx != pdf_memo.end()) {
        return idx->second;
    }

    /* binary search for the nearest item in vals */
    int a = 0;
    int b = n - 1;
    int mid;
    while (a < b) {
        mid = (a + b) / 2;
        if (vals[mid] < x) a = mid + 1;
        else               b = mid;
    }

    static const float eps = 1e-6;
    float ans = 0.0;
    float y;

    /* walk to the left */
    for (int i = a - 1; i >= 0; --i) {
        y = wf_geom(x, vals[i], w) * (float) lens[i];
        ans += y;
        if (y < eps) break;
    }

    /* walk to the right */
    for (int i = a; i < (int) n; ++i) {
        y = wf_geom(x, vals[i], w) * (float) lens[i];
        ans += y;
        if (y < eps) break;
    }

    ans /= (float) m;

    if (!finite(ans)) ans = 0.0;

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

    /* This is the correct way to do this, but it is O(n^2),
     * so we approximate by returning the cdf of the unsmoothed distribution .*/
#if 0

    float ans = 0.0;
    unsigned int u;
    for (u = 0; u <= x; ++u) {
        ans += eval(u);
    }
    return ans;
#endif

    float ans = 0.0;
    size_t i;
    for (i = 0; i < n && vals[i] <= x; ++i) {
        ans += (float) lens[i];
    }
    ans /= (float) m;

    cdf_memo.insert(std::make_pair(x, ans));

    return ans;
}


