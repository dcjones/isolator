
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include "fragment_model.hpp"
#include "sparse_mat.hpp"
#include "transcripts.hpp"

class WeightMatrix;


/* This is the sampler, Isolator's warp-core, so to speak.
 */
class Sampler
{
    public:
        Sampler(const char* bam_fn, const char* ref_fn,
                TranscriptSet& ts, FragmentModel& fm);
        ~Sampler();

        void run(unsigned int num_samples);

    private:
        /* Compute all the entries in frag_probs, given the current value of
         * tmix. */
        void init_frag_probs();

        TranscriptSet& ts;
        FragmentModel& fm;

        /* Transcript mixture coefficients. */
        double* tmix;

        /* Component mixture coefficients. */
        double* cmix;

        WeightMatrix* weight_matrix;
        float* transcript_weights;

        /* transcript_component[i] gives the component number of transcript i */
        unsigned int* transcript_component;

        /* Number of connected components in the fragment/transcript overlap
         * graph. */
        size_t num_components;

        /* Number of transcripts in each component. */
        unsigned int* component_num_transcripts;
        unsigned int** component_transcripts;

        /* frag_counts[i][j] holds the number of occurances of the jth fragment
         * in component i. */
        float** frag_counts;

        /* Row sums of frag_counts, used to compute gradients when optimizing
         * over component mixtures. */
        float* frag_count_sums;

        /* frag_probs[i][j] holds the probability of the jth fragment in
         * component i, given the component and transcript mixtures. */
        float** frag_probs;

        /* A mirror of frag_probs for non-desctructively evaluating proposals.
         * */
        float** frag_probs_prop;

        /* component_frag[i] given the index of the first fragment in component
         * i */
        unsigned int* component_frag;

        friend class InferenceThread;
        friend class MaxPostThread;
        friend class MCMCThread;
        friend double transcript_pair_objf(unsigned int n, const double* xs,
                                           double* grad, void* params);
        friend double component_pair_objf(unsigned int n, const double* xs,
                                          double* grad, void* params);
};


#endif
