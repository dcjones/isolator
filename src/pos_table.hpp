
#ifndef ISOLATOR_POS_TABLE
#define ISOLATOR_POS_TABLE

#include <vector>

#include "common.hpp"
#include "samtools/sam.h"


/* A simple read position structure. */
struct ReadPos
{
    ReadPos(SeqName seqname, uint32_t strand, int32_t pos, uint32_t count)
        : seqname(seqname)
        , strand(strand)
        , pos(pos)
        , count(count)
    {}

    ReadPos()
        : seqname("")
        , strand(strand_na)
        , pos(0)
        , count(0)
    {}

    ReadPos(const ReadPos& other)
        : seqname(other.seqname)
        , strand(other.strand)
        , pos(other.pos)
        , count(other.count)
    {}

    SeqName  seqname;
    uint32_t strand;
    int32_t  pos;
    uint32_t count;
};


/* A class that counts the number of reads starting at a genomic position, used
 * to choose training examples for seqbias.
 *
 * This assumes that reads are added in sorted order. */
class PosTable
{
    public:
        PosTable();
        ~PosTable();

        size_t size();

        /* Tally the start position of the given read. */
        void add(bam1_t* b, samfile_t* bam_f);

        /* Dump to a flat array. */
        void dump(std::vector<ReadPos>& positions, size_t max_size);

    private:
        struct PosTableVal
        {
            pos_t pos;
            unsigned int count;
        };

        struct PosSubtable
        {
            SeqName seqname;
            int32_t tid; /* bam sequence id */

            /* indexed by strand */
            std::vector<PosTableVal> values[2];
        };

        std::vector<PosSubtable> subtables;
};


#endif

