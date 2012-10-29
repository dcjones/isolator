
#ifndef ISOLATOR_TRANSCRIPTS_HPP
#define ISOLATOR_TRANSCRIPTS_HPP

#include <boost/flyweight.hpp>
#include <cstdio>
#include <set>
#include <string>

#include "common.hpp"

/* An exon, relative to some transcript. */
struct Exon
{
    Exon();
    Exon(pos_t, pos_t);

    pos_t start;
    pos_t end;

    /* Order by (start, then end) position. */
    bool operator < (const Exon&) const;
};


/* A transcript is an ordered set of exons with an associated sequence name,
 * strand, etc.
 */
class Transcript : public std::set<Exon>
{
    public:
        Transcript();
        ~Transcript();

        bool operator < (const Transcript& other) const;

        /* Insert a new exon. Use this rather than the std::set insert functions
         * as this will update min_start, max_end. */
        void add(pos_t start, pos_t end);

    private:
        struct gene_id_tag {};
        boost::flyweight<std::string, boost::flyweights::tag<gene_id_tag> > gene_id;

        struct transcript_id_tag {};
        boost::flyweight<std::string, boost::flyweights::tag<transcript_id_tag> > transcript_id;

        struct seq_name_tag {};
        boost::flyweight<std::string, boost::flyweights::tag<seq_name_tag> > seq_name;

        strand_t strand;

        pos_t min_start;
        pos_t max_end;

        pos_t start_codon;
        pos_t stop_codon;

        friend class TranscriptSet;
};


/* A set of transcripts. */
class TranscriptSet
{
    public:
        TranscriptSet();

        /* Read transcripts form a GTF file into the set.
         *
         * Args:
         *   f: A file, opened for reading, containg GTF data.
         */
        void read_gtf(FILE* f);

        /* Number of transcripts held in the set. */
        size_t size() const;

    private:
        /* Transcripts ordered by position. */
        std::set<Transcript> transcripts;
};

#endif


