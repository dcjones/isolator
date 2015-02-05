
#ifndef ISOLATOR_TRANSCRIPTS_HPP
#define ISOLATOR_TRANSCRIPTS_HPP

#include <boost/flyweight.hpp>
#include <cstdio>
#include <set>
#include <string>

#include "common.hpp"
#include "intervals.hpp"
#include "seqbias/twobitseq.hpp"


/* An exon, relative to some transcript. */
struct Exon
{
    Exon();
    Exon(pos_t, pos_t);

    pos_t start;
    pos_t end;

    /* Order by (start, then end) position. */
    bool operator < (const Exon&) const;

    pos_t length() const;
};


/* A transcript is an ordered set of exons with an associated sequence name,
 * strand, etc.
 */
class Transcript : public std::set<Exon>
{
    public:
        Transcript();
        Transcript(const Transcript& other);
        ~Transcript();

        bool operator < (const Transcript& other) const;

        /* Return a reference to the first exon, in terms of positive strand
         * coordinates. */
        Exon& front();

        /* Return a reference to the last exon, in terms of positive strand
         * coordinates. */
        Exon& back();

        bool overlaps(SeqName seqname, pos_t start, pos_t end) const;

        /* Insert a new exon. Use this rather than the std::set insert functions
         * as this will update min_start, max_end. */
        void add(pos_t start, pos_t end);

        /* Return the length of the spliced transcritpt sequence. */
        pos_t exonic_length() const;

        /* Position of the transcript's start site */
        pos_t tss_position() const;

        /* Three prime UTR */
        Interval three_prime_utr() const;

        /* Extract the sequence of the spliced transcript from the full
         * reference sequence.
         *
         * Args:
         *   out: Sequence is written to this object.
         *   ref: The refence sequence the transcript is one. (e.g. chromosome
         *         sequence).
         *   lpad, rpad: Extend the extracted this far to the left and right in
         *               the reference sequence.
         * */
        void get_sequence(twobitseq& out, const twobitseq& ref,
                          pos_t lpad, pos_t rpad) const;

        /* Get the offset of a genomic position within the spliced transcript
         * sequence. */
        pos_t get_offset(pos_t, pos_t leftclip, pos_t rightclip) const;

        GeneID gene_id;
        GeneName gene_name;
        TranscriptID transcript_id;
        SeqName seqname;

        strand_t strand;

        pos_t min_start;
        pos_t max_end;

        pos_t start_codon;
        pos_t stop_codon;

        GeneBiotype biotype;
        GeneSource source;

        /* Transcription group, grouping transcripts whose transcription rate is
         * presumed to be similarly regulated. By default, transcripts with the
         * same transcription start site have the some tgroup number. */
        unsigned int tgroup;

        /* A sequential identifier, unique within the container TranscriptSet. */
        unsigned int id;

    private:
        friend class TranscriptSet;
        friend class TranscriptIntronExonIterator;
};


/* Comparator to order transcripts by transcription start site. */
struct TranscriptCmpTSS
{
    bool operator()(const Transcript&, const Transcript& b) const;
};


/* Comparator to order transcripts by gene_id. */
struct TranscriptCmpGeneId
{
    bool operator()(const Transcript&, const Transcript& b) const;
};


enum IntronExonType {
    EXONIC_INTERVAL_TYPE,
    INTRONIC_INTERVAL_TYPE
};


/* A special iterator that iterates through a transcripts introns and exons in
 * order. */
class TranscriptIntronExonIterator :
    public boost::iterator_facade<TranscriptIntronExonIterator,
                                  const std::pair<Exon, IntronExonType>,
                                  boost::forward_traversal_tag>
{
    public:
        TranscriptIntronExonIterator();
        TranscriptIntronExonIterator(const Transcript&);

    private:
        void increment();
        bool equal(const TranscriptIntronExonIterator&) const;
        const std::pair<Exon, IntronExonType>& dereference() const;

        const Transcript* owner;
        Transcript::const_iterator i;
        std::pair<Exon, IntronExonType> interval;
        friend class boost::iterator_core_access;

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
         *   tgroup_max_tss_dist: distance used to cluster transcripts by TSS
         *                        if use_tss is true.
         *   use_tss: If true define tgroups by tss, otherwise by gene_id.
         */
        void read_gtf(const char* filename, pos_t tgroup_max_tss_dist,
                      bool use_tss,
                      const char* feature="exon",
                      const char* tid_attr="transcript_id",
                      const char* gid_attr="gene_id");

        /* Read annotated alternative exons from a BED file. */
        void read_bed(const char* filename);

        /* Number of transcripts held in the set. */
        size_t size() const;

        /* Return the number of tgroups */
        size_t num_tgroups() const;

        /* Return a vector for each tgroup containing its constituent tids. */
        std::vector<std::vector<unsigned int> > tgroup_tids() const;

        /* Fill a vector with the union of all exonic regions, with strand considered. */
        void get_exonic(std::vector<Interval>&);

        /* Fill a vector with all genomic regions which for every trancsript are
         * either entirely exonic or non-overlapping. */
        void get_consensus_exonic(std::vector<Interval>&);

        /* Fill a vector with intergenic regions, including those at the very
         * beginnings or ends of chromosomes. */
        void get_intergenic(std::vector<Interval>&);

        /* Fill a vector with 5'/3' most exons that don't overlap other
         * transcripts. */
        void get_distinct_5p_3p_exons(const std::vector<Interval>& consensus_exons,
                                      std::vector<Interval>& consensus_5p_exons,
                                      std::vector<Interval>& consensus_3p_exons);

        /* Get 3' UTRs, at lest as long as the given length, and truncated if
         * longer. */
        void get_three_prime_utrs(std::vector<Interval>& intervals, pos_t len);

        /* Find all cassette exons, along with transcripts that exclude or
         * include them. */
        void get_cassette_exons(
                std::vector<Interval>& cassette_exons,
                std::vector<std::vector<unsigned int> >& including_tids,
                std::vector<std::vector<unsigned int> >& excluding_tids);

        void get_retained_introns(
                std::vector<Interval>& retained_introns,
                std::vector<std::vector<unsigned int> >& including_tids,
                std::vector<std::vector<unsigned int> >& excluding_tids);

        /* Acess to std::set::iterator */
        typedef std::set<Transcript>::const_iterator iterator;
        iterator begin();
        iterator end();

    private:
        /* Transcripts ordered by position. */
        std::set<Transcript> transcripts;

        /* Number of transcription groups. */
        size_t _num_tgroups;

        friend class TranscriptSetIterator;
        friend class TranscriptSetLocusIterator;
};


/* A subset of a parent TranscriptSet consisting of all transcripts within a
 * particular locus. */
class TranscriptSetLocus : public std::deque<Transcript>
{
    public:
        TranscriptSetLocus();
        TranscriptSetLocus(const TranscriptSetLocus& other);

        void push_back(Transcript const& t);
        void clear();

        SeqName seqname;
        pos_t min_start;
        pos_t max_end;
};


/* Iterate over overlapping clumps of transcripts.
 *
 * More precisely:
 * Consider an undirected graph where each transcript is a node and there is an
 * edge between two transcripts iff their exonic portions overlap at all. This
 * classs iterates over connected components within that graph.
 */
class TranscriptSetLocusIterator :
    public boost::iterator_facade<TranscriptSetLocusIterator,
                                  const TranscriptSetLocus,
                                  boost::forward_traversal_tag>
{
    public:
        TranscriptSetLocusIterator();
        TranscriptSetLocusIterator(const TranscriptSet&);

    private:
        friend class boost::iterator_core_access;

        void increment();
        bool equal(const TranscriptSetLocusIterator& other) const;
        const TranscriptSetLocus& dereference() const;

        TranscriptSetLocus locus;
        const TranscriptSet* ts;
        std::set<Transcript>::const_iterator i;
};


#endif


