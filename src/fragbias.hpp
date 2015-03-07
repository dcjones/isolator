
#ifndef ISOLATOR_FRAGBIAS
#define ISOLATOR_FRAGBIAS

#include <vector>

#include "common.hpp"

class FragmentModel;

class FragBias
{
    public:
        FragBias(FragmentModel* fm);
        void get_bias(float* out, pos_t l);
        double p;

    private:
        void compute(float* out, pos_t from, pos_t to, pos_t l);

        // Bias is evaluated exactly for a few varying length. Then for specific
        // lengths, we estimate my interpolation.
        std::vector<std::vector<float> > bias;
        std::vector<pos_t> lengths;

        std::vector<float> non_frag_run_pr_sum;
        std::vector<float> non_frag_run_pr;

        pos_t max_non_frag_run;
        FragmentModel* fm;
};


#endif

