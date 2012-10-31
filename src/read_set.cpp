
#include "read_set.hpp"


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


ReadSet::ReadSet()
{
    rs = hattrie_create();
}


ReadSet::~ReadSet()
{
    hattrie_iter_t* i;
    for (i = hattrie_iter_begin(rs, false);
         !hattrie_iter_finished(i);
         hattrie_iter_next(i)) {
        AlignedRead** r = reinterpret_cast<AlignedRead**>(hattrie_iter_val(i));
        delete *r;
    }

    hattrie_iter_free(i);
    hattrie_free(rs);

    std::vector<Alignment*>::iterator j;
    for (j = as.begin(); j != as.end(); ++j) {
        delete *j;
    }
}


void ReadSet::add_alignment(const bam1_t* b)
{
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


