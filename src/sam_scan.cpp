
#include "logger.hpp"
#include "sam_scan.hpp"

extern "C" {
#include "samtools/khash.h"
#include "samtools/sam.h"
KHASH_MAP_INIT_STR(s, int)
}

AlnCountTrie::AlnCountTrie()
{
    t = hattrie_create();
}


AlnCountTrie::~AlnCountTrie()
{
    hattrie_free(t);
}


void AlnCountTrie::inc_mate1(const char* id)
{
    unsigned long* val = hattrie_get(t, id, strlen(id));
    unsigned long cnt = ((*val & 0xffff) + 1) & 0xffff;
    *val = (*val & 0xffff0000) | cnt;
}


void AlnCountTrie::inc_mate2(const char* id)
{
    unsigned long* val = hattrie_get(t, id, strlen(id));
    unsigned long cnt = ((*val >> 16) + 1) & 0xffff;
    *val = (cnt << 16) | (*val & 0xffff);
}


MateCount AlnCountTrie::get(const char* id) const
{
    unsigned long* val = hattrie_tryget(t, id, strlen(id));
    return val == NULL ? MateCount(0, 0) :
                         MateCount(*val & 0xffff, (*val >> 16) & 0xffff);
}


SamScanInterval::SamScanInterval()
    : strand(strand_na)
    , start(-1)
    , end(-1)
    , tid(-1)
{
}


bool SamScanInterval::operator < (const SamScanInterval& other) const
{
    if      (tid != other.tid)     return tid < other.tid;
    else if (start != other.start) return start < other.start;
    else                           return end < other.end;
}


void sam_scan(std::vector<SamScanInterval*>& intervals,
              AlnCountTrie& T, const char* bam_fn, const char* fa_fn)
{
    samfile_t* bam_f;
    bam_f = samopen(bam_fn, "rb", NULL);
    if (bam_f == NULL) {
        bam_f = samopen(bam_fn, "r", NULL);
    }
    if (bam_f == NULL) {
        Logger::abort("Can't open SAM/BAM file [%s].\n", bam_fn);
    }

    /* Sort the intervals in the same order as the BAM file. */
    bam_init_header_hash(bam_f->header);
    khash_t(s)* tbl = reinterpret_cast<khash_t(s)*>(bam_f->header->hash);

    std::vector<SamScanInterval*>::iterator i;
    khiter_t k;
    for (i = intervals.begin(); i != intervals.end(); ++i) {
        k = kh_get(s, tbl, (*i)->seqname.get().c_str());
        if (k == kh_end(tbl)) (*i)->tid = -1;
        else (*i)->tid = kh_value(tbl, k);
    }

    sort(intervals.begin(), intervals.end());

    /* Read the reads. */
    bam1_t* b = bam_init1();
    int32_t last_tid = -1;
    int32_t last_pos = -1;
    while (samread(bam_f, b) > 0) {
        if (b->core.flag & BAM_FUNMAP || b->core.tid < 0) continue;

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort(
                    "Excuse me, but I must insist that your SAM/BAM file be sorted. "
                    "Please run: 'samtools sort'.");
        }
        last_pos = b->core.pos;

        /* Count numbers of alignments by read. */
        if (b->core.flag & BAM_FREAD2) T.inc_mate2(bam1_qname(b));
        else                           T.inc_mate1(bam1_qname(b));

        /* Load sequences as we go. */
        if (b->core.tid != last_tid) {
            /* TODO: load sequence */

            last_tid = b->core.tid;
            last_pos = -1;
        }

        /* Add reads to intervals in which they are contained. */
        /* TODO */
    }

    bam_destroy1(b);
    samclose(bam_f);
}



