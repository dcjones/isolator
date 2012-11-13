
#ifndef ISOLATOR_SAMPLER_HPP
#define ISOLATOR_SAMPLER_HPP

#include "fragment_model.hpp"
#include "sparse_mat.hpp"
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

        /* Fragment weight matrix. */
        SparseMat* W;
};


#endif

