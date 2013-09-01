
#include <boost/foreach.hpp>
#include <deque>

#include "constants.hpp"
#include "gtf/gtf_parse.h"
#include "transcripts.hpp"
#include "trie.hpp"
#include "logger.hpp"


Exon::Exon()
    : start(-1)
    , end(-1)
{
}


Exon::Exon(pos_t start, pos_t end)
    : start(start)
    , end(end)
{
}


bool Exon::operator < (const Exon& other) const
{
    if (this->start == other.start) return this->end < other.end;
    else return this->start < other.start;
}


Transcript::Transcript()
    : strand(strand_na)
    , min_start(-1)
    , max_end(-1)
    , start_codon(-1)
    , stop_codon(-1)
{

}


Transcript::Transcript(const Transcript& other)
    : std::set<Exon>(other)
    , gene_id(other.gene_id)
    , transcript_id(other.transcript_id)
    , seqname(other.seqname)
    , strand(other.strand)
    , min_start(other.min_start)
    , max_end(other.max_end)
    , start_codon(other.start_codon)
    , stop_codon(other.stop_codon)
    , biotype(other.biotype)
    , source(other.source)
    , tgroup(other.tgroup)
    , id(other.id)
{
}


Transcript::~Transcript()
{

}


void Transcript::add(pos_t start, pos_t end)
{
    if (min_start == -1 || min_start > start) min_start = start;
    if (max_end   == -1 || max_end   < end)   max_end   = end;

    std::set<Exon>::insert(Exon(start, end));
}


pos_t Transcript::exonic_length() const
{
    pos_t l = 0;
    for (const_iterator exon = begin(); exon != end(); ++exon) {
        l += exon->end - exon->start + 1;
    }
    return l;
}


void Transcript::get_sequence(twobitseq& dest, const twobitseq& src,
                              pos_t lpad, pos_t rpad) const
{
    const pos_t len = exonic_length();

    dest.resize(len + lpad + rpad);
    pos_t off = 0;

    pos_t tstart = begin()->start;
    while (tstart - lpad < 0) {
        dest.setnuc(off, 0);
        --lpad;
        ++off;
    }

    dest.copy(src, tstart - lpad, off, lpad);
    off += lpad;

    for (const_iterator e = begin(); e != end(); ++e) {
        dest.copy(src, e->start, off, e->end - e->start + 1);
        off += e->end - e->start + 1;
    }

    dest.copy(src, rbegin()->end + 1, off, rpad);
}


pos_t Transcript::get_offset(pos_t pos) const
{
    pos_t offset = 0;
    for (const_iterator exon = begin(); exon != end(); ++exon) {
        if (pos > exon->end) {
            offset += exon->end - exon->start + 1;
        }
        else if (pos >= exon->start) {
            return offset + pos - exon->start;
        }
        // position does not overlap an exon
        else {
            return -1;
        }
    }

    return -1;
}


bool Transcript::operator < (const Transcript& other) const
{
    int c;
    if ((c = seqname_compare(this->seqname.get().c_str(),
                             other.seqname.get().c_str()))) return c < 0;

    else if (min_start != other.min_start) return min_start     < other.min_start;
    else if (max_end   != other.max_end)   return max_end       < other.max_end;
    else if (strand    != other.strand)    return strand        < other.strand;
    else if (gene_id   != other.gene_id)   return gene_id       < other.gene_id;
    else                                   return transcript_id < other.transcript_id;
}


bool Transcript::overlaps(SeqName seqname, pos_t start, pos_t end) const
{
    return this->seqname == seqname &&
           this->min_start <= end &&
           this->max_end   >= start;
}


TranscriptIntronExonIterator::TranscriptIntronExonIterator()
    : owner(NULL)
{
}


TranscriptIntronExonIterator::TranscriptIntronExonIterator(const Transcript& t)
    : owner(&t)
    , i(t.begin())
{
    if (i != owner->end()) {
        interval.first.start = i->start;
        interval.first.end   = i->end;
        interval.second = EXONIC_INTERVAL_TYPE;
    }
}

void TranscriptIntronExonIterator::increment()
{
    if (i == owner->end()) return;

    if (interval.second == EXONIC_INTERVAL_TYPE) {
        pos_t last_end = interval.first.end;
        ++i;
        if (i != owner->end()) {
            interval.first.start = last_end + 1;
            interval.first.end   = i->start - 1;
            interval.second = INTRONIC_INTERVAL_TYPE;
        }
    }
    else {
        interval.first.start = i->start;
        interval.first.end   = i->end;
        interval.second = EXONIC_INTERVAL_TYPE;
    }
}


bool TranscriptIntronExonIterator::equal(const TranscriptIntronExonIterator& other) const
{
    if (owner == NULL || i == owner->end()) {
        return other.owner == NULL || other.i == other.owner->end();
    }
    else if (other.owner == NULL || other.i == other.owner->end()) {
        return false;
    }
    else {
        return i == other.i;
    }
}


const std::pair<Exon, IntronExonType>&TranscriptIntronExonIterator::dereference() const
{
    return interval;
}


/* Construct a TranscriptSet */
TranscriptSet::TranscriptSet()
{
}


size_t TranscriptSet::size() const
{
    return transcripts.size();
}


std::vector<std::vector<unsigned int> > TranscriptSet::tgroup_tids() const
{
    std::vector<std::vector<unsigned int> > tgroup_tids(num_tgroups);
    BOOST_FOREACH (const Transcript& t, transcripts) {
        tgroup_tids[t.tgroup].push_back(t.id);
    }

    return tgroup_tids;
}


void TranscriptSet::read_gtf(FILE* f)
{
    const char* task_name = "Parsing GTF";

    long fsize = 0;
    if (fseek(f, 0, SEEK_END) == 0) {
        fsize = ftell(f);
        rewind(f);
    }

    Logger::push_task(task_name, fsize / 1e6);


    gtf_file_t* gtf_file = gtf_file_alloc(f);
    gtf_row_t* row = gtf_row_alloc();

    size_t count = 0;
    size_t skip_count = 0;

    /* Transcripts indexed by transcript_id. */
    TrieMap<Transcript> ts;

    long fpos_mark = 1e6;

    while (gtf_next(gtf_file, row)) {
        if (fsize > 0 && ftell(f) > fpos_mark) {
            Logger::get_task(task_name).inc();
            fpos_mark += 1e6;
        }

        ++count;

        if (strcmp("exon",        row->feature->s) != 0 &&
            strcmp("start_codon", row->feature->s) != 0 &&
            strcmp("stop_codon",  row->feature->s) != 0) continue;

        str_t* t_id = reinterpret_cast<str_t*>(
                str_map_get(row->attributes, "transcript_id", 13));

        str_t* g_id = reinterpret_cast<str_t*>(
                str_map_get(row->attributes, "gene_id", 7));

        if (t_id == NULL || g_id == NULL) {
            ++skip_count;
            continue;
        }

        str_t* t_biotype = reinterpret_cast<str_t*>(
                str_map_get(row->attributes, "gene_biotype", 12));

        Transcript& t = ts[t_id->s];

        if (t.empty()) {
            t.seqname = row->seqname->s;
            t.gene_id = g_id->s;
            t.transcript_id = t_id->s;
            t.strand = (strand_t) row->strand;
            t.source = row->source->s;
            if (t_biotype) t.biotype = t_biotype->s;
        }

        pos_t pos = (t.strand == strand_pos ? row->start : row->end) - 1;

        if (strcmp("exon", row->feature->s) == 0) {
            /* '-1' to make 0-based, end-inclusive. */
            t.add(row->start - 1, row->end - 1);
        }
        else if (strcmp("start_codon", row->feature->s) == 0) {
            t.start_codon = pos;
        }
        else if (strcmp("stop_codon", row->feature->s) == 0) {
            t.stop_codon = pos;
        }
    }

    gtf_row_free(row);
    gtf_file_free(gtf_file);


    // set transcripts ids and default transcription groups
    typedef std::map<Interval, std::vector<unsigned int> > TSMap;
    typedef std::pair<const Interval, std::vector<unsigned int> > TSMapItem;
    TSMap tss_group;

    std::vector<Transcript> sorted_transcripts;
    sorted_transcripts.reserve(ts.size());

    for (TrieMapIterator<Transcript> t(ts);
         t != TrieMapIterator<Transcript>();
         ++t) {
        sorted_transcripts.push_back(*t->second);
    }
    std::sort(sorted_transcripts.begin(), sorted_transcripts.end());

    unsigned int next_id = 0;
    BOOST_FOREACH (Transcript& t, sorted_transcripts) {
        t.id = next_id++;
        pos_t tss_pos = t.strand == strand_pos ?  t.min_start : t.max_end;
        Interval tss(t.seqname, tss_pos, tss_pos, t.strand);
        tss_group[tss].push_back(t.id);
    }

    unsigned int tgroup = 0;
    BOOST_FOREACH (const TSMapItem& i, tss_group) {
        BOOST_FOREACH (const unsigned int& j, i.second) {
            sorted_transcripts[j].tgroup = tgroup;
        }
        ++tgroup;
    }
    num_tgroups = tgroup;

    // optionally extend the ends of transcripts before inserting them into the
    // set
    BOOST_FOREACH (Transcript& t, sorted_transcripts) {
        /* Extend each transcript 3' and 5' end by some fixed ammount. */
        pos_t u, v;
        Transcript::iterator e1 = t.begin();
        u = e1->start;
        v = e1->end;
        if (t.strand == strand_pos) {
            u = std::max<pos_t>(0, u - constants::transcript_5p_extension);
        }
        else {
            u = std::max<pos_t>(0, u - constants::transcript_3p_extension);
        }
        t.erase(e1);
        t.insert(Exon(u, v));
        t.min_start = u;

        Transcript::reverse_iterator e2 = t.rbegin();
        u = e2->start;
        v = e2->end;
        if (t.strand == strand_pos) {
            v += constants::transcript_3p_extension;
        }
        else {
            v += constants::transcript_5p_extension;
        }
        t.erase(--e2.base());
        t.insert(Exon(u, v));
        t.max_end = v;

        transcripts.insert(t);
    }

    Logger::pop_task(task_name);
    Logger::info("Transcripts: %lu", (unsigned long) size());
}


/* Given a list of exons on one strand, reduce it to a list of disjoint
 * intervals each of which is exonic but overlaps no transcripts ends nor splice
 * junctions. */
static void reduce_exons(std::deque<Exon>& exons)
{
    size_t n = exons.size();
    pos_t* starts = new pos_t[n];
    pos_t* ends   = new pos_t[n];

    size_t j = 0;
    std::deque<Exon>::iterator i;
    for (i = exons.begin(), j = 0; i != exons.end(); ++i, ++j) {
        starts[j] = i->start;
        ends[j]   = i->end;
    }

    std::sort(starts, starts + n);
    std::sort(ends, ends + n);

    std::deque<Exon> c_exons;
    pos_t last = 0;
    int start_count = 0;
    size_t u, v;
    for (u = 0, v = 0; u < n || v < n; ) {

        if (u < n && starts[u] < ends[v]) {
            if (start_count > 0 && last < starts[u]) {
                c_exons.push_back(Exon(last, starts[u] - 1));
            }

            last = std::max(last, starts[u]);
            start_count++;
            u++;
        }
        else if (v < n) {
            if (start_count > 0 && last <= ends[v]) {
                c_exons.push_back(Exon(last, ends[v]));
            }

            last = std::max(last, ends[v] + 1);
            start_count--;
            v++;
        }
    }

    delete [] starts;
    delete [] ends;

    exons = c_exons;
}


void TranscriptSet::get_consensus_exonic(std::vector<Interval>& intervals)
{
    std::set<Transcript>::iterator i;
    Transcript::iterator j;

    /* map strand -> seqname -> exons */
    std::map<strand_t, std::map<SeqName, std::deque<Exon> > > exons;
    std::map<strand_t, std::map<SeqName, std::deque<Exon> > >::iterator k;
    std::map<SeqName, std::deque<Exon> >::iterator l;

    /* construct a list of every exon, indexed by strand and sequence name. */
    for (i = transcripts.begin(); i != transcripts.end(); ++i) {
        k = exons.find((*i).strand);
        if (k == exons.end()) {
            k = exons.insert(
                    std::make_pair(i->strand, std::map<SeqName, std::deque<Exon> >())).first;
        }

        l = k->second.find(i->seqname);
        if (l == k->second.end()) {
            l = k->second.insert(std::make_pair(i->seqname, std::deque<Exon>())).first;
        }

        for (j = i->begin(); j != i->end(); ++j) {
            l->second.push_back(*j);
        }
    }

    /* sort and reduce each list */
    for (k = exons.begin(); k != exons.end(); ++k) {
        for (l = k->second.begin(); l != k->second.end(); ++l) {
            reduce_exons(l->second);
        }
    }

    /* turn the map into a list of intervals. */
    std::deque<Exon>::iterator m;
    for (k = exons.begin(); k != exons.end(); ++k) {
        for (l = k->second.begin(); l != k->second.end(); ++l) {
            for (m = l->second.begin(); m != l->second.end(); ++m) {
                intervals.push_back(Interval(
                                    l->first.get().c_str(),
                                    m->start,
                                    m->end,
                                    k->first));
            }
        }
    }
}


void TranscriptSet::get_exonic(std::vector<Interval>& intervals)
{
    for (strand_t strand = strand_pos; strand != strand_neg; strand = strand_neg) {
        // construct a sorted vector of all exons on the current strand
        std::vector<Interval> exons;
        for (std::set<Transcript>::iterator t = transcripts.begin();
            t != transcripts.end(); ++t) {
            if (t->strand != strand) continue;
            for (Transcript::iterator e = t->begin(); e != t->end(); ++e) {
                exons.push_back(Interval(t->seqname.get().c_str(),
                                         e->start,
                                         e->end,
                                         strand));
            }
        }
        std::sort(exons.begin(), exons.end());

        // find the union
        std::vector<Interval>::iterator i, j;
        for (i = exons.begin(); i != exons.end(); ++i) {
            pos_t start = i->start, end = i->end;
            for (j = i; j != exons.end(); ++j) {
                if (j->start > end) break;
                end = std::max<pos_t>(end, j->end);
            }

            intervals.push_back(Interval(i->seqname.get().c_str(),
                                         start, end, strand));
        }
    }
}


void TranscriptSet::get_intergenic(std::vector<Interval>& intervals)
{
    std::set<Transcript>::iterator i;

    SeqName seqname;
    pos_t end = 0;

    for (i = transcripts.begin(); i != transcripts.end(); ++i) {
        if  (i->seqname != seqname) {
            seqname = i->seqname;
            end = i->max_end;
        }
        else {
            if (i->min_start <= end + 1) {
                end = std::max(end, i->max_end);
            }
            else {
                if (end > 0) {
                    intervals.push_back(
                            Interval(seqname.get().c_str(),
                                     end + 1, i->min_start - 1, strand_na));
                }

                end = i->max_end;
            }
        }
    }
}


void TranscriptSet::get_distinct_5p_3p_exons(const std::vector<Interval>& consensus_exons,
                                             std::vector<Interval>& consensus_5p_exons,
                                             std::vector<Interval>& consensus_3p_exons)
{
    // Group transcripts by gene_id
    std::map<std::string, Transcript> genes;
    std::map<std::string, Transcript>::iterator g;
    std::set<Transcript>::iterator t;
    for (std::set<Transcript>::iterator t = transcripts.begin();
         t != transcripts.end(); ++t) {
        g = genes.find(t->gene_id);
        if (g == genes.end()) {
             g = genes.insert(std::make_pair(t->gene_id.get(), Transcript())).first;
        }

        Transcript::iterator e;
        for (e = t->begin(); e != t->end(); ++e) {
            g->second.insert(*e);
        }
        g->second.strand = t->strand;
        g->second.gene_id = t->gene_id;
        g->second.seqname = t->seqname;
    }

    // Build vectors of 5' and 3' exons.
    std::vector<Interval> exons_5p, exons_3p;

    for (g = genes.begin(); g != genes.end(); ++g) {
        Transcript& t = g->second;

        Transcript::iterator         i_first = t.begin();
        Transcript::reverse_iterator i_last  = t.rbegin();

        Interval first_exon(t.seqname.get().c_str(),
                            i_first->start,
                            i_first->end,
                            t.strand);

        Interval last_exon(t.seqname.get().c_str(),
                           i_last->start,
                           i_last->end,
                           t.strand);

        if (t.strand == strand_pos) {
            exons_5p.push_back(first_exon);
            exons_3p.push_back(last_exon);
        }
        else {
            exons_5p.push_back(last_exon);
            exons_3p.push_back(first_exon);
        }
    }

    std::sort(exons_5p.begin(), exons_5p.end());
    std::sort(exons_3p.begin(), exons_3p.end());

    // find consensus 5' exons
    std::vector<Interval>::const_iterator i = consensus_exons.begin(),
                                          j = exons_5p.begin();

    while (i != consensus_exons.end() && j != exons_5p.end()) {
        while (i != consensus_exons.end() && *i < *j) {
            ++i;
        }
        if (i == consensus_exons.end()) break;

        if (*i == *j) {
            consensus_5p_exons.push_back(*j);
        }
        ++j;
    }

    // find consensus 3' exons
    i = consensus_exons.begin();
    j = exons_3p.begin();

    while (i != consensus_exons.end() && j != exons_3p.end()) {
        while (i != consensus_exons.end() && *i < *j) {
            ++i;
        }
        if (i == consensus_exons.end()) break;

        if (*i == *j) {
            consensus_3p_exons.push_back(*j);
        }
        ++j;
    }
}


TranscriptSet::iterator TranscriptSet::begin()
{
    return transcripts.begin();
}


TranscriptSet::iterator TranscriptSet::end()
{
    return transcripts.end();
}


TranscriptSetLocus::TranscriptSetLocus()
    : min_start(-1)
    , max_end(-1)
{
}


TranscriptSetLocus::TranscriptSetLocus(const TranscriptSetLocus& other)
    : std::deque<Transcript>(other)
    , seqname(other.seqname)
    , min_start(other.min_start)
    , max_end(other.max_end)
{
}


void TranscriptSetLocus::push_back(Transcript const& t)
{
    if (min_start == -1 || t.min_start < min_start) min_start = t.min_start;
    if (max_end   == -1 || t.max_end   > max_end)   max_end   = t.max_end;
    seqname = t.seqname;
    std::deque<Transcript>::push_back(t);
}


void TranscriptSetLocus::clear()
{
    min_start = -1;
    max_end = -1;
    std::deque<Transcript>::clear();
}


TranscriptSetLocusIterator::TranscriptSetLocusIterator()
    : ts(NULL)
{
}


TranscriptSetLocusIterator::TranscriptSetLocusIterator(const TranscriptSet& ts)
    : ts(&ts)
{
    i = ts.transcripts.begin();
    increment();
}


void TranscriptSetLocusIterator::increment()
{
    if (ts == NULL) return;

    locus.clear();
    for (; i != ts->transcripts.end(); ++i) {
        if (locus.empty() ||
                i->overlaps(locus.seqname, locus.min_start, locus.max_end)) {
            locus.push_back(*i);
        }
        else break;
    }

    /* Mark the end of the iterator. */
    if (locus.empty()) ts = NULL;
}


bool TranscriptSetLocusIterator::equal(const TranscriptSetLocusIterator& other) const
{
    return (ts == NULL && other.ts == NULL) || i == other.i;
}


const TranscriptSetLocus& TranscriptSetLocusIterator::dereference() const
{
    return locus;
}


