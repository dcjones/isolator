
#include "read_set.hpp"
#include "constants.hpp"


Alignment::Alignment()
    : start(-1)
    , end(-1)
    , cigar_len(0)
    , cigar(NULL)
    , strand(strand_na)
{
}


Alignment::Alignment(const Alignment& a)
{
    start     = a.start;
    end       = a.end;
    strand    = a.strand;
    cigar_len = a.cigar_len;
    cigar = new uint32_t[cigar_len];
    memcpy(cigar, a.cigar, cigar_len * sizeof(uint32_t));
}


Alignment::Alignment(const bam1_t* b)
{
    start     = (pos_t)b->core.pos;
    end       = bam_calend(&b->core, bam1_cigar(b)) - 1;
    strand    = bam1_strand(b);
    cigar_len = b->core.n_cigar;
    cigar     = new uint32_t [b->core.n_cigar];
    memcpy(cigar, bam1_cigar(b), cigar_len * sizeof(uint32_t));
}


Alignment::~Alignment()
{
    delete cigar;
}


bool Alignment::operator == (const bam1_t* b) const
{
    if (this->start != (pos_t) b->core.pos) return false;
    if (this->strand != bam1_strand(b)) return false;
    if (this->cigar_len != b->core.n_cigar) return false;
    if (this->end != (pos_t) (bam_calend(&b->core, bam1_cigar(b)) - 1)) return false;
    return memcmp(cigar, bam1_cigar(b), cigar_len * sizeof(uint32_t)) == 0;
}


bool Alignment::operator != (const bam1_t* b) const
{
    return !(*this == b);
}


AlignedRead::AlignedRead()
{
}


AlignedRead::~AlignedRead()
{
}


bool AlignedRead::operator < (const AlignedRead& other) const
{
    if (this->mate1.size() != other.mate1.size()) {
        return this->mate1.size() < other.mate1.size();
    }

    if (this->mate2.size() != other.mate2.size()) {
        return this->mate2.size() < other.mate2.size();
    }

    if      (start != other.start) return start < other.start;
    else if (end   != other.end)   return end   < other.end;


    std::set<Alignment*> S, T;

    S.insert(this->mate1.begin(), this->mate1.end());
    T.insert(other.mate1.begin(), other.mate1.end());

    if (S != T) return S < T;

    S.clear();
    T.clear();

    S.insert(this->mate2.begin(), this->mate2.end());
    T.insert(other.mate2.begin(), other.mate2.end());

    return S < T;
}


bool AlignmentPair::valid_frag() const
{
    /* Reads with only mate1 are considered valid single-end reads, for our
     * purposes */
    if (mate1 == NULL) return false;
    if (mate1 != NULL && mate1 == NULL) return true;

    switch (constants::libtype) {
        case constants::LIBTYPE_FR:
            if (mate1 ->strand == mate2->strand) return false;
            if (mate1->strand == strand_pos)     return mate1->start <= mate2->start;
            else                                 return mate1->start >= mate2->start;

        case constants::LIBTYPE_RF:
            if (mate1->strand == mate2->strand)  return false;
            if (mate1->strand == strand_pos)     return mate1->start >= mate2->start;
            else                                 return mate1->start <= mate2->start;

        case constants::LIBTYPE_FF:
            if (mate1->strand != mate2->strand) return false;
            if (mate1->strand == strand_pos)    return mate1->start <= mate2->start;
            else                                return mate1->start >= mate2->start;

        default:
            return false;
    }
}


pos_t AlignmentPair::naive_frag_len() const
{
    if (mate1 == NULL || mate2 == NULL) return 0;
    return std::max(mate2->end - mate1->start, mate1->end - mate2->start);
}


AlignedReadIterator::AlignedReadIterator()
    : r(NULL)
    , i(0)
    , j(0)
{
}


AlignedReadIterator::AlignedReadIterator(const AlignedRead& r)
    : r(&r)
    , i(0)
    , j(0)
{
    if (i < r.mate1.size()) p.mate1 = r.mate1[i];
    if (j < r.mate2.size()) p.mate2 = r.mate2[j];
}


AlignedReadIterator::~AlignedReadIterator()
{
}


void AlignedReadIterator::increment()
{
    if (finished()) return;
    do {
        if (r->paired && j < r->mate2.size()) {
            do {
                ++j;
            } while (j < r->mate2.size() && r->mate2[j] == r->mate2[j - 1]);
        }

        if (!r->paired || j >= r->mate2.size()) {
            do {
                ++i;
            } while(i < r->mate1.size() && r->mate1[i] == r->mate1[i - 1]);

            j = 0;
        }

        p.mate1 = i < r->mate1.size() ? r->mate1[i] : NULL;
        p.mate2 = (j < r->mate2.size() && r->paired) ? r->mate2[j] : NULL;
    } while (!finished() && !p.valid_frag());
}


bool AlignedReadIterator::finished() const
{
    return r == NULL || i >= r->mate1.size() || (r->paired && j >= r->mate2.size());
}


bool AlignedReadIterator::equal(const AlignedReadIterator& other) const
{
    if (finished()) {
        return other.finished();
    }
    else {
        if (other.finished()) return false;
        return r == other.r && i == other.i && j == other.j;
    }
}


const AlignmentPair& AlignedReadIterator::dereference() const
{
    return p;
}


ReadSet::ReadSet()
    : rs(NULL)
{
}


ReadSet::~ReadSet()
{
    clear();

    std::vector<Alignment*>::iterator j;
    for (j = as.begin(); j != as.end(); ++j) {
        delete *j;
    }
}


void ReadSet::add_alignment(const bam1_t* b)
{
    if (rs == NULL) rs = hattrie_create();

    AlignedRead** rp = reinterpret_cast<AlignedRead**>(
            hattrie_get(rs, bam1_qname(b), b->core.l_qname - 1));

    if (*rp == NULL) *rp = new AlignedRead();
    AlignedRead* r = *rp;

    if (as.empty() || *as.back() != b) as.push_back(new Alignment(b));
    Alignment* a = as.back();

    if (b->core.flag & BAM_FREAD1 || b->core.flag & BAM_FREAD2) r->paired = true;

    if (b->core.flag & BAM_FREAD2) r->mate2.push_back(a);
    else                           r->mate1.push_back(a);

    if (r->start == -1 || r->start < a->start) r->start = a->start;
    if (r->end   == -1 || r->end   < a->end)   r->end    = a->end;
}


void ReadSet::clear()
{
    if (rs == NULL) return;

    hattrie_iter_t* i;
    for (i = hattrie_iter_begin(rs, false);
         !hattrie_iter_finished(i);
         hattrie_iter_next(i)) {
        AlignedRead** r = reinterpret_cast<AlignedRead**>(hattrie_iter_val(i));
        delete *r;
    }

    hattrie_iter_free(i);
    hattrie_free(rs);
    rs = NULL;
}


void ReadSet::make_unique_read_counts(ReadSet::UniqueReadCounts& counts)
{
    if (rs == NULL) return;
    hattrie_iter_t* i;
    ReadSet::UniqueReadCounts::iterator j;
    for (i = hattrie_iter_begin(rs, false);
         !hattrie_iter_finished(i);
         hattrie_iter_next(i)) {
        AlignedRead** r = reinterpret_cast<AlignedRead**>(hattrie_iter_val(i));

        j = counts.find(*r);
        if (j == counts.end()) {
            counts.insert(std::make_pair(*r, 1));
        }
        else j->second += 1;
    }
}



