
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include "fragment_model.hpp"
#include "queue.hpp"
#include "sample_db.hpp"
#include "sparse_mat.hpp"
#include "transcripts.hpp"
#include <vector>

class WeightMatrix;
class AbundanceSamplerThread;
class MultireadSamplerThread;


// Unit of work used in the abundance sampler.
struct ComponentBlock
{
    ComponentBlock()
        : u(INT_MAX)
        , v(INT_MAX)
    {
    }

    ComponentBlock(unsigned int u, unsigned int v)
        : u(u) , v(v)
    {
    }

    bool is_end_of_queue() const
    {
        return u == INT_MAX;
    }

    unsigned int u, v;
};


// Unit of work used in the multiread sampler.
struct MultireadBlock
{
    MultireadBlock()
        : u(INT_MAX)
        , v(INT_MAX)
    {
    }

    MultireadBlock(unsigned int u, unsigned int v)
        : u(u), v(v)
    {
    }

    bool is_end_of_queue() const
    {
        return u == INT_MAX;
    }

    unsigned int u, v;
};


/* This is the sampler, Isolator's warp-core, so to speak.
 */
class Sampler
{
    public:
        Sampler(const char* bam_fn, const char* ref_fn,
                TranscriptSet& ts, FragmentModel& fm);
        ~Sampler();

        void run(unsigned int num_samples, SampleDB& out, bool run_gc_correction);

    private:
        // Run a single multiread sampler round.
        void sample_multireads();

        // Run a single tmix/cmix round.
        void sample_abundance();

        /* Compute all the entries in frag_probs, given the current value of
         * tmix. */
        void init_frag_probs();

        // Zero the count of each multiread alignment
        void init_multireads();

        // Recompute frag_count_sums
        void update_frag_count_sums();

        /* Perform post-hoc adjustment for transcript GC-content. */
        void gc_correction(float* xs, size_t num_samples);

        TranscriptSet& ts;
        FragmentModel& fm;

        // Multiread sampling
        std::vector<MultireadSamplerThread*> multiread_threads;
        Queue<MultireadBlock> multiread_queue;

        // Abundance sampling
        std::vector<AbundanceSamplerThread*> abundance_threads;
        Queue<ComponentBlock> component_queue;

        /* Transcript mixture coefficients. */
        double* tmix;

        /* Component mixture coefficients. */
        double* cmix;

        WeightMatrix* weight_matrix;
        float* transcript_weights;

        /* GC content of each transcript. */
        float* transcript_gc;

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

        /* Number of multireads. */
        unsigned num_multireads;

        /* Number of alignments for each multiread. */
        unsigned int* multiread_num_alignments;

        /* A single multiread alignment. */
        struct MultireadAlignment
        {
            unsigned int component;
            unsigned int frag;
            float align_pr;
        };

        /* A flat array of all the multiread alignments */
        MultireadAlignment* multiread_alignment_pool;

        /* Pointers into multiread_alignment_poll for each multiread. */
        MultireadAlignment** multiread_alignments;

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

        /* Samples indexed by transcript id */
        float** samples;

        friend class InferenceThread;
        friend class AbundanceSamplerThread;
        friend class MultireadSamplerThread;
};


#endif
