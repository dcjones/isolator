
#ifndef ISOLATOR_READ_SET_HPP
#define ISOLATOR_READ_SET_HPP

#include <vector>

#include "common.hpp"
#include "hat-trie/hat-trie.h"
#include "samtools/sam.h"

/* A representation of an aligned sequence. */
struct Alignment
{
    Alignment();
    Alignment(const Alignment& other);
    Alignment(const bam1_t* other);
    ~Alignment();

    bool operator == (const bam1_t* b) const;
    bool operator != (const bam1_t* b) const;

    pos_t start;
    pos_t end;
    uint16_t  cigar_len;
    uint32_t* cigar;
    uint8_t   strand;
};


/* A read with some number of alignments. */
struct AlignedRead
{
    AlignedRead();
    ~AlignedRead();

    pos_t start, end;
    bool paired;
    std::vector<Alignment*> mate1;
    std::vector<Alignment*> mate2;
};


/* A container for a set of reads indexed by id. */
class ReadSet
{
    public:
        ReadSet();
        ~ReadSet();

        /* Add an alignment to the read set. */
        void add_alignment(const bam1_t* b);

    private:
        /* Map of read ids to a AligneRead objects. */
        hattrie_t* rs;

        /* A pool of unique alignments. */
        std::vector<Alignment*> as;
};

#endif

