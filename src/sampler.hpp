
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include <vector>

#include "fragment_model.hpp"
#include "gc.hpp"
#include "queue.hpp"
#include "sample_db.hpp"
#include "sparse_mat.hpp"
#include "transcripts.hpp"

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
                TranscriptSet& ts, FragmentModel& fm,
                bool run_gc_correction, bool use_priors=false);
        ~Sampler();

        void run(unsigned int num_samples, SampleDB& out);

        // Called prior to any calls to iterate
        void start();

        // Generate a new sample. start() must be called prior to running an
        // iterations.
        void sample();

        // Return the current sampler state: a vector containing the relative
        // abundance of each transcript indexed by tid.
        const std::vector<double>& state() const;

        // Scaling factor to convert between relative abundance, normalized by sequencing depth,
        // estimates in terms of number of reads.
        const std::vector<double>& expr_scaling() const;

        // Use priors on transcript start site usage and splicing
        void engage_priors();
        void disengage_priors();

        // Hyperparameters that this sampler is not responsible for updating.
        // Public so they can be updated without a hassle.
        struct {
            // Normalization factor accounting for sequencing depth
            double scale;

            // Degrees-of-freedom parameter for tgroup abundance prior
            double tgroup_nu;

            // Location parameter for tgroup abundance prior
            std::vector<double> tgroup_mu;

            // Scale parameter for the tgroup abundance prior
            std::vector<double> tgroup_sigma;

            // TODO: Hyperparameters for splicing. Yikes.
        } hp;

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

        TranscriptSet& ts;
        FragmentModel& fm;
        GCCorrection* gc_correct;

        // Multiread sampling
        std::vector<MultireadSamplerThread*> multiread_threads;
        Queue<MultireadBlock> multiread_queue;

        // Abundance sampling
        std::vector<AbundanceSamplerThread*> abundance_threads;
        Queue<ComponentBlock> component_queue;

        // True if priors on cmix and tmix are used. If false, hyperparameters
        // (values in the hp structure) are ignored.
        bool use_priors;

        /* Transcript mixture coefficients. Within-component relative abundance. */
        double* tmix;

        /* Transcription group within-component relative abundance. */
        double* tgroupmix;

        /* Component mixture coefficients. */
        double* cmix;

        WeightMatrix* weight_matrix;
        double* transcript_weights;

        /* Factor by which transcript estimates in terms of number of reads differ
           from the estimate accounting for transcript length, gc bias, and sequencing depth */
        std::vector<double> tgroup_scaling;

        /* GC content of each transcript. */
        double* transcript_gc;

        /* transcript_component[i] gives the component number of transcript i */
        unsigned int* transcript_component;

        /* transcript_tgroups[i] gives the tgroup of transcript i */
        unsigned int* transcript_tgroup;

        /* map tgroup i to a vector of its constituent tids. */
        std::vector<std::vector<unsigned int> > tgroup_tids;

        /* map component i to a vector of its constituent tgroups */
        std::vector<std::set<unsigned int> > component_tgroups;

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

        /* Total number of fragments beign considered  */
        double total_frag_count;

        /* frag_probs[i][j] holds the probability of the jth fragment in
         * component i, given the component and transcript mixtures. */
        float** frag_probs;

        /* A mirror of frag_probs for non-desctructively evaluating proposals.
         * */
        float** frag_probs_prop;

        /* component_frag[i] given the index of the first fragment in component
         * i */
        unsigned int* component_frag;

        /* Temporary space to do normalization on expression before copying to
         * samples */
        std::vector<double> expr;

        /* Samples indexed by transcript id */
        float** samples;

        friend class InferenceThread;
        friend class AbundanceSamplerThread;
        friend class MultireadSamplerThread;
};


#endif
