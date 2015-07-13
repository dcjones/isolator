
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include <vector>

#include "fragment_model.hpp"
#include "queue.hpp"
#include "transcripts.hpp"

class WeightMatrix;
class AbundanceSamplerThread;


// Unit of work used in the abundance sampler.
struct ComponentBlock
{
    ComponentBlock()
        : u(INT_MAX)
        , v(INT_MAX)
        , rng(NULL)
    {
    }

    ComponentBlock(unsigned int u, unsigned int v, rng_t& rng)
        : u(u) , v(v), rng(&rng)
    {
    }

    bool is_end_of_queue() const
    {
        return u == INT_MAX;
    }

    unsigned int u, v;
    rng_t* rng;
};



/* This is the sampler, Isolator's warp-core, so to speak.
 */
class Sampler
{
    public:
        Sampler(unsigned int rng_seed,
                const char* bam_fn, const char* ref_fn,
                std::set<std::string> excluded_seqs,
                TranscriptSet& ts, FragmentModel& fm,
                bool run_frag_correction,
                bool use_priors=false);
        ~Sampler();

        // Called prior to any calls to iterate
        void start();

        // Called when sampling should finish.
        void stop();

        // Generate a new sample. start() must be called prior to running an
        // iterations.
        void sample();

        // Instead drawing a random sample, hillclimp towards the maximum
        // posterior.
        void optimize();

        // Return the current sampler state: a vector containing the relative
        // abundance of each transcript indexed by tid.
        const std::vector<float>& state() const;

        // Use priors on transcript start site usage and splicing
        void engage_priors();
        void disengage_priors();

        // Hyperparameters that this sampler is not responsible for updating.
        // Public so they can be updated without a hassle.
        struct {
            // Normalization factor accounting for sequencing depth
            double scale;

            // Mean parameter for transcript abundance prior
            std::vector<float> mean;

            // Shape parameter for transcript abundance prior
            std::vector<float> shape;

            // Normal splice mean
            std::vector<double> splice_mu;

            // Normal splice sigma
            std::vector<double> splice_sigma;
        } hp;

        unsigned long num_frags() const;

        unsigned long num_alignments() const;

    private:
        // Run a single tmix/cmix round.
        void sample_abundance();

        /* Compute all the entries in frag_probs, given the current value of
         * tmix. */
        void init_frag_probs();

        // Abundance sampling
        std::vector<AbundanceSamplerThread*> abundance_threads;
        Queue<ComponentBlock> component_queue;
        Queue<int> component_notify_queue;
        std::vector<rng_t> abundance_rng_pool;

        // Components are processed in this order to help prevent
        // thread starvation.
        std::vector<unsigned int> ordered_components;

        // True if priors on cmix and tmix are used. If false, hyperparameters
        // (values in the hp structure) are ignored.
        bool use_priors;

        /* Transcript mixture coefficients. Within-component relative abundance. */
        double* tmix;

        /* Component mixture coefficients. */
        double* cmix;

        WeightMatrix* weight_matrix;
        double* transcript_weights;

        /* Factor by which transcript estimates in terms of number of reads differ
           from the estimate accounting for transcript length, gc bias, and sequencing depth */
        std::vector<float> expr_scaling;

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

        /* frag_probs[i][j] holds the probability of the jth fragment in
         * component i, given the component and transcript mixtures. */
        float** frag_probs;

        /* component_frag[i] given the index of the first fragment in component
         * i */
        unsigned int* component_frag;

        /* Temporary space to do normalization on expression before copying to
         * samples */
        std::vector<float> expr;

        friend class AbundanceSamplerThread;
        friend class ComponentSampler;
        friend class InterTranscriptSampler;
};


#endif
