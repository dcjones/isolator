
#include "gtf/gtf_parse.h"
#include "transcripts.hpp"
#include "trie.hpp"


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
    if ((c = seqname_compare(this->seq_name.get().c_str(),
                             other.seq_name.get().c_str()))) return c < 0;

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


void TranscriptSet::read_gtf(FILE* f)
{
    gtf_file_t* gtf_file = gtf_file_alloc(f);
    gtf_row_t* row = gtf_row_alloc();

    size_t count = 0;
    size_t skip_count = 0;

    /* Transcripts indexed by transcript_id. */
    TrieMap<Transcript> ts;

    while (gtf_next(gtf_file, row)) {
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
            t.seq_name = row->seqname->s;
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
       transcripts.insert(t->second); 
    }

}


