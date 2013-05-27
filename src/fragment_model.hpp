

#ifndef ISOLATOR_FRAGMENT_MODEL_HPP
#define ISOLATOR_FRAGMENT_MODEL_HPP

#include "hat-trie/hat-trie.h"
#include "seqbias/sequencing_bias.hpp"
#include "transcripts.hpp"
#include "emp_dist.hpp"

typedef std::pair<unsigned int, unsigned int> MateCount;


/* This is a trie that maps read ids to a pair of integers giving the number
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

        /* Set the number of alignments for the read with teh given ID. */
        void set(const char*, const MateCount&);

        /* True if the key is present. */
        bool has(const char*) const;

        /* Number of entries. */
        size_t size() const;

    private:
        hattrie_t* t;

        friend class AlnCountTrieIterator;
};


class AlnCountTrieIterator :
    public boost::iterator_facade<AlnCountTrieIterator,
                                  const std::pair<const char*, MateCount>,
                                  boost::forward_traversal_tag>
{
    public:
        AlnCountTrieIterator();
        AlnCountTrieIterator(const AlnCountTrie&);
        ~AlnCountTrieIterator();

    private:
        friend class boost::iterator_core_access;
        void increment();
        bool equal(const AlnCountTrieIterator& other) const;
        const std::pair<const char*, MateCount>& dereference() const;

        hattrie_iter_t* it;

        /* key value at the current position. */
        std::pair<const char*, MateCount> x;
};


/* Assign indexes to a set of string. (Read ids in this case.) */
class AlnIndex
{
    public:
        AlnIndex();
        ~AlnIndex();

        size_t size() const;
        void clear();

        void add(const char* key);

        /* Return -1 if the key is not present, otherwise return the key's
         * index. */
        int get(const char* key);


    private:
        hattrie_t* t;
        boost::mutex mut;
};


/* A probabalistic model of fragment sampling in RNA-Seq experiments. */
class FragmentModel
{
    public:
        FragmentModel();
        ~FragmentModel();
        void estimate(TranscriptSet& ts, const char* bam_fn, const char* fa_fn);
        void train_seqbias(TranscriptSet& ts, const char* bam_fn, const char* fa_fn);

        /* Fragment length probability ,cdf, and median using a fallback when no
         * emperical distribution is available. */
        float frag_len_p(pos_t frag_len);
        float frag_len_c(pos_t frag_len);
        float frag_len_med();

        /* On it's pass through the reads, estimate will also count the number
         * of alignments for each read. */
        AlnIndex multireads;

        /* A set of reads that have a very large number of alignments, that if
         * included would dramatically increase the run time. */
        AlnIndex blacklist;

        /* Probability of a read being from the same strand as the transcript it
         * originates from. */
        float strand_specificity;

        /* A model of sequence bias. */
        sequencing_bias* sb[3];

        /* Distribution over fragment lengths. */
        EmpDist* frag_len_dist;

        /* Three-prime/five-prime distribution indexed by length bin and strand.
         * */
        std::vector<double> tp_dist[4][2];
};


#endif

