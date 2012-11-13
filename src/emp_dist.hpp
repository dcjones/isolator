
#ifndef ISOLATOR_EMP_DIST
#define ISOLATOR_EMP_DIST

#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <climits>


/* A smoothed emperical distribution over [0, inf] using parzen windows. */
class EmpDist
{
    public:
        EmpDist(EmpDist&);

        /* Construct a emperical distribution from n observations stored in xs.
         *
         * Input in run-length encoded samples, which must be in sorted order.
         *
         * Args:
         *   vals: An array of unique observations.
         *   lens: Nuber of occurances for each observation.
         *   n: Length of vals and lens.
         *   w: Smoothing parameter. A larger number corresponds to more
         *      smoothing.
         */
        EmpDist(const unsigned int* vals,
                const unsigned int* lens,
                size_t n,
                float w = 0.1);

        ~EmpDist();

        /* Compute the median of the distribution. */
        float median() const;

        /* Probability density function. */
        float pdf(unsigned int x) const;

        /* Comulative p robabability. */
        float cdf(unsigned int x) const;

    private:
        /* Unprotected pdf computation.  */
        float eval(unsigned int x) const;

        /* Run-length encoded samples. (That is, a histogram.) */
        unsigned int* vals;
        unsigned int* lens;

        /* Length of the `vals` and `lens` arrays. */
        size_t n;

        /* Sum of `lens`. */
        size_t m;

        /* Smoothing coefficient */
        float w;

        /* Precomputed median */
        mutable float med;

        /* Memoization of the pdf and cdf functions. */
        mutable boost::unordered_map<unsigned int, float> pdf_memo;
        mutable boost::unordered_map<unsigned int, float> cdf_memo;

        /* Mutex for pdf_memo/cdf_memo access. */
        mutable boost::mutex mut;
};



#endif


