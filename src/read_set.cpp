
#include <cctype>
#include "read_set.hpp"
#include "constants.hpp"


Alignment::Alignment()
    : start(-1)
    , end(-1)
    , cigar_len(0)
    , cigar(NULL)
    , strand(strand_na)
    , alnnum(0)
{
}


Alignment::Alignment(const Alignment& a)
{
    start     = a.start;
    end       = a.end;
    strand    = a.strand;
    alnnum    = a.alnnum;
    misaligned_pr  = a.misaligned_pr;
    cigar_len = a.cigar_len;
    cigar = new uint32_t[cigar_len];
    memcpy(cigar, a.cigar, cigar_len * sizeof(uint32_t));
}


Alignment::Alignment(const bam1_t* b)
{
    start     = bam_truepos(&b->core, bam1_cigar(b));
    end       = bam_trueend(&b->core, bam1_cigar(b)) - 1;
    strand    = bam1_strand(b);

    unsigned int mismatch_count = 0;
    char* s = reinterpret_cast<char*>(bam_aux_get(b, "NM"));
    if (!s) s = reinterpret_cast<char*>(bam_aux_get(b, "nM"));
    if (s) {
        mismatch_count = (unsigned int) bam_aux2i(reinterpret_cast<unsigned char*>(s));
    }

    s = reinterpret_cast<char*>(bam_aux_get(b, "MD"));
    double misalign_likelihood[2] = {1.0, 1.0};
    if (s) {
        unsigned char* qs = bam1_qual(b);
        pos_t off = 0;
        char* t = s;
        unsigned long skip;
        while (*s) {
            skip = strtoul(s, &t, 10);
            if (t > s) {
                off += (pos_t) skip;
            }

            // probability the base is miscalled
            double p = pow(10.0, - (double) qs[off] / 10);
            misalign_likelihood[0] *=
                p + (1.0 - p) * constants::misalign_mismatch_pr;
            misalign_likelihood[1] *=
                p + (1.0 - p) * constants::mismatch_pr;

            if (!*t) break;

            s = t + 1;
            ++off;
        }

    }
    else {
        misalign_likelihood[0] = pow(constants::misalign_mismatch_pr, mismatch_count);
        misalign_likelihood[1] = pow(constants::mismatch_pr, mismatch_count);
    }

        misaligned_pr = (constants::misalign_prior * misalign_likelihood[0]) /
            ((constants::misalign_prior * misalign_likelihood[0]) +
             (1.0 - constants::misalign_prior) * misalign_likelihood[1]);


    s = reinterpret_cast<char*>(bam_aux_get(b, "HI"));
    alnnum = s ? (uint16_t) bam_aux2i(reinterpret_cast<unsigned char*>(s)) : 0;

    cigar_len = b->core.n_cigar;
    cigar     = new uint32_t [b->core.n_cigar];
    memcpy(cigar, bam1_cigar(b), cigar_len * sizeof(uint32_t));
}


Alignment::~Alignment()
{
    delete [] cigar;
}


bool Alignment::operator == (const bam1_t* b) const
{
    if (this->start != (pos_t) bam_truepos(&b->core, bam1_cigar(b))) return false;
    if (this->strand != bam1_strand(b)) return false;

    uint16_t alnnum = 0;
    unsigned char* s = bam_aux_get(b, "HI");
    if (s) alnnum = (uint16_t) bam_aux2i(s);
    if (this->alnnum != alnnum) return false;

    if (this->cigar_len != b->core.n_cigar) return false;
    if (this->end != (pos_t) (bam_trueend(&b->core, bam1_cigar(b)) - 1)) return false;
    return memcmp(cigar, bam1_cigar(b), cigar_len * sizeof(uint32_t)) == 0;
}


bool Alignment::operator != (const bam1_t* b) const
{
    return !(*this == b);
}


bool Alignment::operator < (const Alignment& other) const
{
    if      (start  != other.start)  return start < other.start;
    else if (end    != other.end)    return end < other.end;
    else if (strand != other.strand) return strand < other.strand;
    else if (alnnum != other.alnnum) return alnnum < other.alnnum;
    //else if (mapq   != other.mapq)   return mapq < other.mapq;
    //else if (mismatch_count != other.mismatch_count) return mismatch_count < other.mismatch_count;
    else if (cigar_len != other.cigar_len) return cigar_len < other.cigar_len;
    else {
        return memcmp(cigar, other.cigar, cigar_len * sizeof(uint32_t)) < 0;
    }
}


bool Alignment::operator == (const Alignment& other) const
{
    return start     == other.start &&
           end       == other.end &&
           strand    == other.strand &&
           alnnum    == other.alnnum &&
           //mapq      == other.mapq &&
           cigar_len == other.cigar_len &&
           memcmp(cigar, other.cigar, cigar_len * sizeof(uint32_t)) == 0;
}


CigarIterator::CigarIterator()
    : owner(NULL)
    , i(0)
{
}


CigarIterator::CigarIterator(const Alignment& owner)
    : owner(&owner)
    , i(0)
{
    if (i < owner.cigar_len) {
        c.start = owner.start;
        c.end   = c.start + (owner.cigar[i] >> BAM_CIGAR_SHIFT) - 1;
        c.op    = owner.cigar[i] & BAM_CIGAR_MASK;
    }
    else {
        c.start = -1;
        c.end   = -1;
        c.op    = BAM_CPAD;
    }
}


void CigarIterator::increment()
{
    if (i >= owner->cigar_len) return;

    if (++i < owner->cigar_len) {
        c.start = c.end + 1;
        c.end   = c.start + (owner->cigar[i] >> BAM_CIGAR_SHIFT) - 1;
        c.op    = owner->cigar[i] & BAM_CIGAR_MASK;
    }
}


bool CigarIterator::equal(const CigarIterator& other) const
{
    if (owner == NULL || i >= owner->cigar_len) {
        return other.owner == NULL || other.i >= other.owner->cigar_len;
    }
    else if (other.owner == NULL || other.i >= other.owner->cigar_len) {
        return false;
    }
    else {
        return owner == other.owner && i == other.i;
    }
}


const Cigar& CigarIterator::dereference() const
{
    return c;
}


AlignedRead::AlignedRead()
    : start(-1)
    , end(-1)
    , paired(false)
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

    std::set<Alignment*, AlignmentPtrCmp> S, T;
    std::set<Alignment*, AlignmentPtrCmp>::iterator i, j;

    if (this->mate1.size() > 0) {
        S.insert(this->mate1.begin(), this->mate1.end());
        T.insert(other.mate1.begin(), other.mate1.end());

        if (S.size() != T.size()) return S.size() < T.size();

        i = S.begin();
        j = T.begin();
        while (i != S.end()) {
            if (**i < **j) return true;
            ++i;
            ++j;
        }

        S.clear();
        T.clear();
    }

    if (this->mate2.size() > 0) {
        S.insert(this->mate2.begin(), this->mate2.end());
        T.insert(other.mate2.begin(), other.mate2.end());

        if (S.size() != T.size()) return S.size() < T.size();

        i = S.begin();
        j = T.begin();
        while (i != S.end()) {
            if (**i < **j) return true;
            ++i;
            ++j;
        }
    }

    return false;
}


AlignmentPair::AlignmentPair()
    : mate1(NULL)
    , mate2(NULL)
{
}


AlignmentPair::AlignmentPair(const AlignmentPair& other)
    : mate1(other.mate1)
    , mate2(other.mate2)
{
}


bool AlignmentPair::operator == (const AlignmentPair& other) const
{
    if (mate1 == other.mate2 ||
        (mate1 != NULL && other.mate1 != NULL && *mate1 == *other.mate1)) {

        return mate2 == other.mate2 ||
            (mate2 != NULL && other.mate2 != NULL && *mate2 == *other.mate2);

    }
    else return false;
}


bool AlignmentPair::operator < (const AlignmentPair& other) const
{
    if (mate1 == other.mate1 ||
        (mate1 != NULL && other.mate1 != NULL && *mate1 == *other.mate1)) {

        if (mate2 == NULL)            return other.mate2 != NULL;
        else if (other.mate2 == NULL) return false;
        else                          return *mate2 < *other.mate2;
    }
    else {
        if (mate1 == NULL)            return true;
        else if (other.mate1 == NULL) return false;
        else                          return *mate1 < *other.mate1;
    }
}


bool AlignmentPair::valid_frag() const
{
    if (mate1 == NULL || mate2 == NULL) return true;
    if (mate1->alnnum != mate2->alnnum) return false;

    switch (constants::libtype) {
        case constants::LIBTYPE_FR:
            if (mate1->strand == mate2->strand) return false;
            if (mate1->strand == strand_pos)    return mate1->start <= mate2->start;
            else                                return mate1->start >= mate2->start;

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


std::pair<pos_t, pos_t> AlignmentPair::flanks_soft_clipping() const
{
    std::pair<pos_t, pos_t> clipping(0, 0);

    if (mate1 && mate2 && mate1->cigar_len > 0 && mate2->cigar_len > 0) {
        const Alignment *a1, *a2;
        if (mate1->strand == strand_pos) {
            a1 = mate1;
            a2 = mate2;
        }
        else {
            a1 = mate2;
            a2 = mate1;
        }

        if ((a1->cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
            clipping.first = a1->cigar[0] >> BAM_CIGAR_SHIFT;
        }

        if ((a2->cigar[a2->cigar_len - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
            clipping.second = a2->cigar[a2->cigar_len - 1] >> BAM_CIGAR_SHIFT;
        }
    }
    else {
        const Alignment* a = mate1 ? mate1 : mate2;
        if (a->cigar_len > 0) {
            if (a->strand == strand_pos &&
                    (a->cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
                clipping.first = a->cigar[0] >> BAM_CIGAR_SHIFT;
            }
            else if ((a->cigar[a->cigar_len - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
                clipping.second = a->cigar[a->cigar_len - 1] >> BAM_CIGAR_SHIFT;
            }
        }
    }

    return clipping;
}


/* A couple functions to assist with AlignmentPair::frag_len */
static bool exon_compatible_cigar_op(uint8_t op)
{
    return op == BAM_CMATCH || op == BAM_CSOFT_CLIP || op == BAM_CINS || op == BAM_CDEL;
}


static bool intron_compatible_cigar_op(uint8_t op)
{
    return op == BAM_CREF_SKIP || op == BAM_CSOFT_CLIP;
}


pos_t AlignmentPair::frag_len(const Transcript& t) const
{
    /* Reorder mates for convenience. */
    const Alignment *a1, *a2;
    if (mate1 != NULL && mate2 != NULL) {
        if (mate1->start <= mate2->start) {
            a1 = mate1;
            a2 = mate2;
        }
        else {
            a1 = mate2;
            a2 = mate1;
        }

#if 0
        /* NOTE: we currently don't allow one mate to be contained within the
         * other. This might occur if the mates are if different lengths, but
         * isn't likely to be a problem. This function will need some
         * modifications if it does prove to be an issue. */
        if (a1->end > a2->end) {
            return -1;
        }
#endif
    }
    else if (mate1 != NULL) {
        a1 = mate1;
        a2 = NULL;
    }
    else {
        a1 = mate2;
        a2 = NULL;
    }

    if (a1 && a1->start < t.min_start) {
        return -1;
    }
    if (a2 && a2->end > t.max_end) {
        return -1;
    }

    TranscriptIntronExonIterator e1(t);
    CigarIterator c1(*a1);
    pos_t intron_len = 0;

    // Skip any leading soft clip
    if (c1 != CigarIterator() && c1->op == BAM_CSOFT_CLIP) ++c1;
    Cigar c = *c1;

    while (e1 != TranscriptIntronExonIterator() && c1 != CigarIterator()) {
        // case 1: e entirely preceedes c
        if (c.start > e1->first.end) {
            ++e1;
        }

        // case 2: c is contained within e
        else if (c.end   >= e1->first.start
              && c.end   <= e1->first.end
              && c.start >= e1->first.start)
        {
            if (e1->second == EXONIC_INTERVAL_TYPE) {
                if (!exon_compatible_cigar_op(c.op)) {
                    return -1;
                }
            }
            else {
                if (!intron_compatible_cigar_op(c.op)) {
                    return -1;
                }
                intron_len += e1->first.end - e1->first.start + 1;
                ++e1;
            }

            ++c1;
            c = *c1;
        }

        // case 3: soft clipping partially overlapping an exon or intron
        else if (c.op == BAM_CSOFT_CLIP) {
            ++c1;
            c = *c1;
        }

        else if (c.end > e1->first.end && c.op == BAM_CMATCH) {
            if (e1->second == EXONIC_INTERVAL_TYPE &&
                c.end - e1->first.end <= constants::max_intron_overlap) {
                c.end = e1->first.end;
            }
            else if (e1->second == INTRONIC_INTERVAL_TYPE &&
                     e1->first.end >= c.start &&
                     e1->first.end - c.start < constants::max_intron_overlap) {
                c.start = e1->first.end + 1;
            }
            else return -1;
        }

        // case 4: c precedes partially overlaps e
        else {
            return -1;
        }
    }

    if (c1 != CigarIterator()) {
        return -1;
    }

    /* alignment is compatible, but single ended. */
    if (a2 == NULL) {
        return 0;
    }

    bool e2_sup_e1 = false; /* marks when e2 > e1 */
    TranscriptIntronExonIterator e2(t);
    CigarIterator c2(*a2);
    c = *c2;

    while (e2 != TranscriptIntronExonIterator() && c2 != CigarIterator()) {
        // case 1: e entirely preceedes c
        if (c.start > e2->first.end) {
            if (e2->second == INTRONIC_INTERVAL_TYPE && e2_sup_e1) {
                intron_len += e2->first.end - e2->first.start + 1;
            }

            if (e1 == e2) e2_sup_e1 = true;

            ++e2;
        }

        // case 2: c is contained within e
        else if (c.end   >= e2->first.start
              && c.end   <= e2->first.end
              && c.start >= e2->first.start)
        {
            if (e2->second == EXONIC_INTERVAL_TYPE) {
                if (!exon_compatible_cigar_op(c.op)) {
                    return -1;
                }
            }
            else {
                if (!intron_compatible_cigar_op(c.op)) {
                    return -1;
                }
            }

            ++c2;
            c = *c2;
        }

        // case 3: soft clipping partially overlapping an exon or intron
        else if (c.op == BAM_CSOFT_CLIP) {
            ++c2;
            c = *c2;
        }

        else if (c.end > e2->first.end && c.op == BAM_CMATCH) {
            if (e2->second == EXONIC_INTERVAL_TYPE &&
                c.end - e2->first.end <= constants::max_intron_overlap) {
                c.end = e2->first.end;
            }
            else if (e2->second == INTRONIC_INTERVAL_TYPE &&
                     e2->first.end >= c.start &&
                     e2->first.end - c.start < constants::max_intron_overlap) {
                c.start = e2->first.end + 1;
            }
            else return -1;
        }

        // case 4: c precedes or partially overlaps e
        else {
            return -1;
        }
    }

    // skip any trailing soft clip
    if (c2 != CigarIterator() && c2->op == BAM_CSOFT_CLIP) ++c2;

    if (c2 != CigarIterator()) {
        return -1;
    }

    pos_t fraglen = a2->end - a1->start + 1 - intron_len;
    assert(fraglen > 0);
    return fraglen;
}


pos_t AlignmentPair::naive_frag_len() const
{
    if (mate1 == NULL || mate2 == NULL) return 0;
    pos_t len = 1 + std::max(mate2->end - mate1->start, mate1->end - mate2->start);
    /* The fragment length must be at leas as long as the reads themselves. */
    if (len < mate2->end - mate2->start + 1 ||
        len < mate1->end - mate1->start + 1) {
        return 0;
    }
    else return len;
}


strand_t AlignmentPair::strand() const
{
    if (mate1) return (strand_t) mate1->strand;
    else {
        switch (constants::libtype) {
            case constants::LIBTYPE_FF:
                return (strand_t) mate2->strand;

            case constants::LIBTYPE_FR:
            case constants::LIBTYPE_RF:
            default:
                return other_strand((strand_t) mate2->strand);
        }
    }
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
    p.mate1 = r.mate1.empty() ? NULL : r.mate1[0];
    p.mate2 = r.mate2.empty() ? NULL : r.mate2[0];
}


AlignedReadIterator::~AlignedReadIterator()
{
}


void AlignedReadIterator::increment()
{
    if (finished()) return;
    do {
        if (r->mate1.size() == 0) {
            do {
                ++j;
            } while (j < r->mate2.size() && r->mate2[j] == r->mate2[j - 1]);
        }
        else if (r->mate2.size() == 0) {
            do {
                ++i;
            } while (i < r->mate1.size() && r->mate1[i] == r->mate1[i - 1]);
        }
        else {
            do {
                ++j;
            } while (j < r->mate2.size() && r->mate2[j] == r->mate2[j - 1]);

            if (j >= r->mate2.size()) {
                j = 0;
                do {
                    ++i;
                } while (i < r->mate1.size() && r->mate1[i] == r->mate1[i - 1]);
            }
        }

        p.mate1 = i < r->mate1.size() ? r->mate1[i] : NULL;
        p.mate2 = j < r->mate2.size() ? r->mate2[j] : NULL;
    } while (!finished() && !p.valid_frag());
}


bool AlignedReadIterator::finished() const
{
    return r == NULL ||
        (r->mate1.size() == 0 && j >= r->mate2.size()) ||
        (r->mate2.size() == 0 && i >= r->mate1.size()) ||
        (r->mate1.size() > 0 && r->mate2.size() > 0 && i >= r->mate1.size());
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

    if (r->start == -1 || r->start > a->start) r->start = a->start;
    if (r->end   == -1 || r->end   < a->end)   r->end   = a->end;
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

    hattrie_iter_free(i);
}


ReadSetIterator::ReadSetIterator()
    : it(NULL)
{
}


ReadSetIterator::ReadSetIterator(const ReadSet& s)
    : it(NULL)
{
    if (s.rs) {
        it = hattrie_iter_begin(s.rs, false);
        if (!hattrie_iter_finished(it)) {
            x.first = hattrie_iter_key(it, NULL);
            x.second = *reinterpret_cast<AlignedRead**>(hattrie_iter_val(it));
        }
    }
}


ReadSetIterator::~ReadSetIterator()
{
    hattrie_iter_free(it);
}


void ReadSetIterator::increment()
{
    hattrie_iter_next(it);
    if (!hattrie_iter_finished(it)) {
        x.first = hattrie_iter_key(it, NULL);
        x.second = *reinterpret_cast<AlignedRead**>(hattrie_iter_val(it));
    }
}


bool ReadSetIterator::equal(const ReadSetIterator& other) const
{
    if (it == NULL || hattrie_iter_finished(it)) {
        return other.it == NULL || hattrie_iter_finished(other.it);
    }
    else if (other.it == NULL || hattrie_iter_finished(other.it)) {
        return false;
    }
    else return hattrie_iter_equal(it, other.it);
}


const std::pair<const char*, AlignedRead*>& ReadSetIterator::dereference() const
{
    return x;
}


