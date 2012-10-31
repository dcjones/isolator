
#ifndef ISOLATOR_SAM_SCAN_HPP
#define ISOLATOR_SAM_SCAN_HPP

#include <vector>

#include "common.hpp"
#include "hat-trie/hat-trie.h"
#include "read_set.hpp"
#include "samtools/sam.h"
#include "intervals.hpp"

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


/* A genomic interval with an associated function.
 *
 * Given an instance of SamScanInterval, sam_scan will index the reads
 * overlapping the interval, and when there are no more reads within the
 * interval, call the finish function.
 *
 * Subclasss should define the finish function, naturally.
 * */
class SamScanInterval
{
    public:
        SamScanInterval();
        SamScanInterval(const Interval& interval);
        void add_alignment(const bam1_t*);
        virtual void finish() = 0;

        /* Comparison based on genomic position. */
        bool operator < (const SamScanInterval& other) const;

        /* Clear used memory. */
        void clear();

    private:
        SeqName seqname;
        pos_t start, end;
        strand_t strand;

        ReadSet rs;

        /* A the sequence ID assigned by the BAM file, so we can arrange
         * intervals in the same order. */
        int32_t tid;

        friend void sam_scan(std::vector<SamScanInterval*>& intervals,
                             AlnCountTrie& T, const char* bam_fn);
};


/* This function makes a single pass over a SAM/BAM file collecting summary
 * statistics.
 *
 * To avoid seeking around in the data or making multiple passes, this system is
 * a bit complicated.
 *
 * Args:
 *   T: Count mate alignments for every read in the given trie.
 *   bam_fn: Filename of a SAM/BAM file.
 * */
void sam_scan(std::vector<SamScanInterval*>& intervals,
              AlnCountTrie& T, const char* bam_fn);

#endif

