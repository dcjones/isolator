
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

    private:
        TranscriptSet& ts;
        FragmentModel& fm;

        WeightMatrix* weight_matrix;
        float* transcript_weights;

        /* transcript_component[i] gives the component number of transcript i */
        unsigned int* transcript_component;

        /* Number of connected components in the fragment/transcript overlap
         * graph. */
        size_t num_components;

        /* frag_couns[i][j] holds the number of occurances of the jth fragment
         * in component i. */
        float** frag_counts;

        /* component_frag[i] given the index of the first fragment in component
         * i */
        unsigned int* component_frag;
};


#endif
