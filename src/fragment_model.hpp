

#ifndef ISOLATOR_FRAGMENT_MODEL_HPP
#define ISOLATOR_FRAGMENT_MODEL_HPP

#include <string>

#include "hat-trie/hat-trie.h"
#include "seqbias/sequencing_bias.hpp"
#include "transcripts.hpp"
#include "emp_dist.hpp"
#include "gcbias.hpp"
#include "tpbias.hpp"


/* Assign indexes to a set of string. (Read ids in this case.) */
class AlnIndex
{
    public:
        AlnIndex();
        ~AlnIndex();

        size_t size() const;
        void clear();

        void add(const char* key, uint16_t alnnum);

        /* Return -1 if the key is not present, otherwise return the key's
         * index. */
        long get(const char* key, uint16_t alnnum);


    private:
        size_t makekey(const char* key, uint16_t alnnum);

        char* keybuf;
        size_t keybuf_size;

        hattrie_t* t;
        boost::mutex mut;
};


/* A probabalistic model of fragment sampling in RNA-Seq experiments. */
class FragmentModel
{
    public:
        FragmentModel();
        ~FragmentModel();
        void estimate(TranscriptSet& ts, const char* bam_fn, const char* fa_fn,
                      bool use_gc_correction, bool use_3p_correction,
                      bool tabulate_seqbias, std::set<std::string> bias_training_seqnames);

        /* Fragment length probability ,cdf, and median using a fallback when no
         * emperical distribution is available. */
        float frag_len_p(pos_t frag_len);
        float frag_len_c(pos_t frag_len);
        float frag_len_med();

        /* Numerical indices assigned to fragment ids */
        AlnIndex alnindex;

        /* Probability of a read being from the same strand as the transcript it
         * originates from. */
        float strand_specificity;

        /* A model of sequence bias. */
        sequencing_bias* sb[2];

        /* Tabulations of sequence bias. */
        SeqbiasTabulation sb_tabulation[2];

        /* A model of fragment GC bias. */
        GCBias* gcbias;

        /* A model of 3'-5' positional bias */
        TPBias* tpbias;

        /* Distribution over fragment lengths. */
        EmpDist* frag_len_dist;

};


#endif

