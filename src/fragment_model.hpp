

#ifndef ISOLATOR_FRAGMENT_MODEL_HPP
#define ISOLATOR_FRAGMENT_MODEL_HPP

#include "hat-trie/hat-trie.h"
#include "seqbias/sequencing_bias.hpp"
#include "transcripts.hpp"
#include "emp_dist.hpp"

typedef std::pair<unsigned int, unsigned int> MateCount;


/* This a is a trie that maps read ids to a pair of integers giving the number
 * of alignments of the first and second mate, respectively. */
class AlnCountTrie
{
    public:
        AlnCountTrie();
        ~AlnCountTrie();

        /* Record an occurance of a mate1 alignment of the given ID. */
        void inc_mate1(const char*);

        /* Record an occurance of a mate2 alignment of the given ID. */
        void inc_mate2(const char*);

        /* Get the number of alignments of both mates of the given ID. */
        MateCount get(const char*) const;

    private:
        hattrie_t* t;
};


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

