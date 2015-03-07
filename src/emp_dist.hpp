
#ifndef ISOLATOR_EMP_DIST
#define ISOLATOR_EMP_DIST

#include <climits>
#include <vector>


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
         */
        EmpDist(const unsigned int* vals,
                const unsigned int* lens,
                size_t n);

        /* Compute the median of the distribution. */
        float median() const;

        /* Probability density function. */
        float pdf(unsigned int x) const;

        /* Comulative p robabability. */
        float cdf(unsigned int x) const;

    private:
        std::vector<float> pdfvals;
        std::vector<float> cdfvals;

        /* Precomputed median */
        float med;
};



#endif


