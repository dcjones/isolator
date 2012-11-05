
#include <deque>

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


Transcript::~Transcript()
{

}


void Transcript::add(pos_t start, pos_t end)
{
    if (min_start == -1 || min_start > start) min_start = start;
    if (max_end   == -1 || max_end   < end)   max_end   = end;

    std::set<Exon>::insert(Exon(start, end));
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



/* Construct a TranscriptSet */
TranscriptSet::TranscriptSet()
{
}


size_t TranscriptSet::size() const
{
    return transcripts.size();
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

        Transcript& t = ts[t_id->s];

        if (t.empty()) {
            t.seqname = row->seqname->s;
            t.gene_id = g_id->s;
            t.transcript_id = t_id->s;
            t.strand = (strand_t) row->strand;
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

    for (TrieMapIterator<Transcript> t(ts);
         t != TrieMapIterator<Transcript>();
         ++t) {
        transcripts.insert(*t->second);
    }

    Logger::pop_task(task_name);
    Logger::info("Read %lu transcripts.", (unsigned long) size());
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


void TranscriptSet::get_intergenic(std::vector<Interval>& intervals)
{
    std::set<Transcript>::iterator i;

    SeqName seqname;
    pos_t start = 0, end = 0;

    for (i = transcripts.begin(); i != transcripts.end(); ++i) {
        if  (i->seqname != seqname) {
            seqname = i->seqname;
            start = 0;
            end = i->max_end;
        }
        else {
            if (i->min_start < end) {
                end = std::max(end, i->min_start);
            }
            else {
                if (end > 0) {
                    intervals.push_back(
                            Interval(seqname.get().c_str(),
                                     end + 1, i->min_start - 1, strand_na));
                }

                start = i->min_start;
                end = i->max_end;
            }
        }
    }
}

