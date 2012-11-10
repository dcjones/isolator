
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include "fragment_model.hpp"
#include "transcripts.hpp"

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

        /* What do we need?
         *
         * A dense vector of mixture coefficients for components.
         *
         * For each component
         *   1. a dense vector of transcript mixture coefficients.
         *   2. a sparse matrix of fragment probabilities, given transcript of
         *      origin.
         *
         *
         * For each component, we need a sparse matrix
         *
         *
         * Paralellization:
         *
         *  It's not strictly necessary to sample components in parallel. We
         *  could run two or more chains seperately.
         *
         *  We would need more memory though. Recall that we are destructively
         *  updating the matrix giving fragment probabilities.
         *
         */

        /* Before we write the actual sampler though, we need to figure out if
         * there is a possible prior on transcript abundance that would allow
         * efficient slicing.
         *
         * What would a power law distribution look like? Can we easily find the
         * gradient?
         * */
};


#endif

