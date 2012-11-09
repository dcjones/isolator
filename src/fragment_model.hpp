

#ifndef ISOLATOR_FRAGMENT_MODEL_HPP
#define ISOLATOR_FRAGMENT_MODEL_HPP

#include "sam_scan.hpp"
#include "seqbias/sequencing_bias.hpp"
#include "transcripts.hpp"
#include "emp_dist.hpp"

/* A probabalistic model of fragment sampling in RNA-Seq experiments. */
class FragmentModel
{
    public:
        FragmentModel();
        ~FragmentModel();
        void estimate(TranscriptSet& ts, const char* bam_fn, const char* fa_fn);

    private:
        /* On it's pass through the reads, estimate will also count the number
         * of alignments for each read. */
        AlnCountTrie alncnt;

        /* Probability of a read being from the same strand as the transcript it
         * originates from. */
        float strand_specificity;

        /* A model of sequence bias. */
        sequencing_bias* sb;

        /* Distribution over fragment lengths. */
        EmpDist* frag_len_dist;

};


#endif

