
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include "constants.hpp"
#include "hat-trie/hat-trie.h"
#include "logger.hpp"
#include "queue.hpp"
#include "read_set.hpp"
#include "sampler.hpp"
#include "samtools/sam.h"
#include "samtools/samtools_extra.h"

extern "C" {
#include "samtools/khash.h"
KHASH_MAP_INIT_STR(s, int)
}


/* A trivial class to serve unique indices to multile threads. */
class Indexer
{
    public:
        Indexer()
            : next(0)
        {
        }

        unsigned int get()
        {
            boost::lock_guard<boost::mutex> lock(mut);
            return next++;
        }

        unsigned int count()
        {
            return next;
        }

    private:
        boost::mutex mut;
        unsigned int next;
};


/* An entry in the coordinate respresentation of the sparse weight matrix. */
struct WeightMatrixEntry
{
    WeightMatrixEntry(unsigned int i, unsigned int j, float w)
        : i(i), j(j), w(w)
    {}

    bool operator < (const WeightMatrixEntry& other) const
    {
        return i != other.i ? i < other.i : j < other.j;
    }

    struct ColOrd
    {
        bool operator () (const WeightMatrixEntry& a, const WeightMatrixEntry& b) const
        {
            return a.j < b.j;
        }
    };

    unsigned int i, j;
    float w;

};


/* Store information needed to process multireads. */
struct Multiread
{
    Multiread()
    {
    }

    /* Information necessary to compute one entry in the likelihood matrix. */
    struct Cell
    {
        Cell(unsigned int transcript_idx, float frag_weight, float align_pr)
            : transcript_idx(transcript_idx)
            , frag_weight(frag_weight)
            , align_pr(align_pr)
        {
        }

        unsigned int transcript_idx;
        float frag_weight;
        float align_pr;
    };

    struct CellTranOrd {
        bool operator () (const Cell& a, const Cell& b) const
        {
            return a.transcript_idx < b.transcript_idx;
        }
    };

    std::deque<Cell> cells;
};


/* When initializing the likelihood matrix rows for reads with multiple
 * alignments cannot be added until all alignments have been read.
 * This class collects information for multireads to be processed later. */
class MultireadIndex :
    public boost::unordered_map<unsigned int, Multiread>
{
    public:
        void push(unsigned int key, const Multiread::Cell& cell)
        {
            Multiread& r = (*this)[key];
            r.cells.push_back(cell);
        }

        void push(unsigned int key, const std::deque<Multiread::Cell>& cells)
        {
            Multiread& r = (*this)[key];
            //r.cells.reserve(r.cells.size() + cells.size());
            r.cells.insert(r.cells.end(), cells.begin(), cells.end());
        }
};


/* An interval of overlapping transcripts, which constitutes a unit of work when
 * doing initialization across multiple threads. */
class SamplerInitInterval
{
    public:
        SamplerInitInterval(
                TranscriptSetLocus ts,
                Queue<SamplerInitInterval*>& q)
            : ts(ts)
            , tid(-1)
            , q(q)
        {
        }

        void add_alignment(const bam1_t* b)
        {
            rs.add_alignment(b);
        }

        void clear()
        {
            rs.clear();
        }

        void finish()
        {
            q.push(this);
        }

        bool operator < (const SamplerInitInterval& other) const
        {
            if (tid != other.tid) return tid < other.tid;
            else if (ts.seqname != other.ts.seqname) {
                return ts.seqname < other.ts.seqname;
            }
            else if (ts.min_start != other.ts.min_start) {
                return ts.min_start < other.ts.min_start;
            }
            else {
                return ts.max_end < other.ts.max_end;
            }
        }

        TranscriptSetLocus ts;
        ReadSet rs;
        boost::shared_ptr<twobitseq> seq;

        int32_t tid;
    private:
        Queue<SamplerInitInterval*>& q;

        friend void sam_scan(std::vector<SamplerInitInterval*>& intervals,
                             const char* bam_fn,
                             const char* fa_fn,
                             const char* task_name);
};


struct SamplerInitIntervalPtrCmp
{
    bool operator () (const SamplerInitInterval* a,
                      const SamplerInitInterval* b)
    {
        return *a < *b;
    }
};


/* Read through a sorted SAM/BAM file, initializing the sampler as we go. */
void sam_scan(std::vector<SamplerInitInterval*>& intervals,
              const char* bam_fn, const char* fa_fn,
              const char* task_name)
{
    /* Measure file size to monitor progress. */
    size_t input_size = 0;
    size_t input_block_size = 1000000;
    if (task_name) {
        FILE* f = fopen(bam_fn, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            input_size = (size_t) ftell(f);
            fclose(f);
        }
        Logger::push_task(task_name, input_size / input_block_size);
    }

    /* Open the SAM/BAM file */
    samfile_t* bam_f;
    bam_f = samopen(bam_fn, "rb", NULL);
    if (bam_f == NULL) {
        bam_f = samopen(bam_fn, "r", NULL);
    }
    if (bam_f == NULL) {
        Logger::abort("Can't open SAM/BAM file %s.", bam_fn);
    }

    /* Open the FASTA file */
    faidx_t* fa_f = NULL;
    if (fa_fn) {
        fa_f = fai_load(fa_fn);
        if (fa_f == NULL) {
            Logger::abort("Can't open FASTA file %s.", fa_fn);
        }
    }

    /* Sort the intervals in the same order as the BAM file. */
    bam_init_header_hash(bam_f->header);
    khash_t(s)* tbl = reinterpret_cast<khash_t(s)*>(bam_f->header->hash);

    std::vector<SamplerInitInterval*>::iterator i;
    khiter_t k;
    for (i = intervals.begin(); i != intervals.end(); ++i) {
        k = kh_get(s, tbl, (*i)->ts.seqname.get().c_str());
        if (k == kh_end(tbl)) (*i)->tid = -1;
        else (*i)->tid = kh_value(tbl, k);
    }

    sort(intervals.begin(), intervals.end(), SamplerInitIntervalPtrCmp());

    /* First interval which the current read may be contained in. */
    size_t j, j0 = 0;
    size_t n = intervals.size();

    size_t last_file_pos = 0, file_pos;
    size_t read_num = 0;

    /* Read the reads. */
    bam1_t* b = bam_init1();
    int32_t last_tid = -1;
    int32_t last_pos = -1;
    while (samread(bam_f, b) > 0) {
        ++read_num;

        if (read_num % 1000 == 0) {
            file_pos = samtell(bam_f);
            if (file_pos >= last_file_pos + input_block_size && input_size > 0) {
                Logger::get_task(task_name).inc();
                last_file_pos = file_pos;
            }
        }

        if (b->core.flag & BAM_FUNMAP || b->core.tid < 0) continue;

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort(
                    "Excuse me, but I must insist that your SAM/BAM file be sorted. "
                    "Please run: 'samtools sort'.");
        }

        if (fa_f && b->core.tid != last_tid) {
            /* TODO: handle the case in which an alignments sequence is not in
             * the reference We should ignore these reads and complain loudly.
             * .*/
            int seqlen;
            char* seqstr = faidx_fetch_seq(
                    fa_f, bam_f->header->target_name[b->core.tid],
                    0, INT_MAX, &seqlen);

            if (seqstr == NULL) {
                Logger::abort("Couldn't read sequence %s",
                        bam_f->header->target_name[b->core.tid]);
            }

            boost::shared_ptr<twobitseq> seq(new twobitseq(seqstr));
            free(seqstr);

            for (j = j0; j < n && b->core.tid == intervals[j]->tid; ++j) {
                intervals[j]->seq = seq;
            }

            last_pos = -1;
        }

        last_tid = b->core.tid;
        last_pos = b->core.pos;

        /* Add reads to intervals in which they are contained. */
        for (j = j0; j < n; ++j) {
            if (b->core.tid < intervals[j]->tid) break;
            if (b->core.tid > intervals[j]->tid) {
                assert(j == j0);
                intervals[j0++]->finish();
                continue;
            }

            if (b->core.pos < intervals[j]->ts.min_start) break;
            if (b->core.pos > intervals[j]->ts.max_end) {
                if (j == j0) {
                    intervals[j0++]->finish();
                }
                continue;
            }

            pos_t b_end = (pos_t) bam_calend(&b->core, bam1_cigar(b)) - 1;
            if (b_end <= intervals[j]->ts.max_end) {
                intervals[j]->add_alignment(b);
            }
        }
    }

    for(; j0 < n; ++j0) intervals[j0]->finish();

    bam_destroy1(b);
    samclose(bam_f);
    if (fa_f) fai_destroy(fa_f);

    if (task_name) Logger::pop_task(task_name);
}


class SamplerInitThread
{
    public:
        SamplerInitThread(
                FragmentModel& fm,
                Indexer& read_indexer,
                Queue<SamplerInitInterval*>& q)
            : fm(fm)
            , read_indexer(read_indexer)
            , q(q)
            , thread(NULL)
        {

        }

        ~SamplerInitThread()
        {
        }

        void run()
        {
            SamplerInitInterval* locus;
            while (true) {
                if ((locus = q.pop()) == NULL) break;
                process_locus(locus);
                delete locus;
            }
        };

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&SamplerInitThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        void process_locus(SamplerInitInterval* locus);

        /* Compute sequence bias for both mates on both strand.
         *
         * Results are stored in mate1_seqbias, mate2_seqbias. */
        void transcript_sequence_bias(const SamplerInitInterval& locus,
                                      const Transcript& t);

        /* Compute the sum of the weights of all fragmens in the transcript.
         *
         * This assumes that transcript_sequence_bias has already been called,
         * and the mate1_seqbias/mate2_seqbias arrays set for t. */
        float transcript_weight(const Transcript& t);

        /* Compute (a number proportional to) the probability of observing the
         * given fragment from the given transcript, assuming every transcript
         * is equally expressed. */
        float fragment_weight(
                const Transcript& t,
                const AlignmentPair& a);

        std::vector<WeightMatrixEntry> weight_matrix_entries;

        typedef std::pair<unsigned int, unsigned int> FragIdxCount;
        std::vector<FragIdxCount> frag_counts;

        typedef std::pair<unsigned int, float> TransIdxCount;
        std::vector<TransIdxCount> transcript_weights;

        FragmentModel& fm;
        MultireadIndex multiread_index;
        Indexer& read_indexer;
        Queue<SamplerInitInterval*>& q;
        boost::thread* thread;

        /* Temprorary space for computing sequence bias, indexed by strand. */
        std::vector<float> mate1_seqbias[2];
        std::vector<float> mate2_seqbias[2];

        /* Temporary space for holding transcript sequences. */
        twobitseq tseq0; /* + strand */
        twobitseq tseq1; /* reverse complement of tseq0 */

        /* Temporary space used when computing transcript weight. */
        std::vector<float> ws;

        friend class Sampler;
};


void SamplerInitThread::process_locus(SamplerInitInterval* locus)
{
    /* A map of fragments onto sequential indices. */
    typedef std::map<AlignmentPair, FragIdxCount> Frags;
    Frags frag_idx;

    /* The set of fragments within this locus that have no compatible
     * transcripts. */
    typedef std::set<AlignmentPair> FragSet;
    FragSet excluded_frags;

    /* The set of reads with multiple alignment. */
    typedef std::set<std::pair<unsigned int, AlignedRead*> > MultireadSet;
    MultireadSet multiread_set;

    /* Collapse identical reads, filter out those that don't overlap any
     * transcript. */
    size_t num_entries = 0;
    for (ReadSetIterator r(locus->rs); r != ReadSetIterator(); ++r) {
        /* Skip multireads for now. */
        int multiread_num = fm.multireads.get(r->first);
        if (multiread_num >= 0) {
            multiread_set.insert(std::make_pair((unsigned int) multiread_num,
                                                r->second));
            continue;
        }

        AlignedReadIterator a(*r->second);
        if (a == AlignedReadIterator()) continue;

        if (excluded_frags.find(*a) != excluded_frags.end()) continue;

        Frags::iterator f = frag_idx.find(*a);
        if (f != frag_idx.end()) {
            f->second.second++;
            continue;
        }

        bool has_compatible_transcript = false;
        for (TranscriptSetLocus::iterator t = locus->ts.begin();
                t != locus->ts.end(); ++t) {
            if (a->frag_len(*t) >= 0) {
                has_compatible_transcript = true;
            }
            ++num_entries;
        }

        if (has_compatible_transcript) {
            FragIdxCount& f = frag_idx[*a];
            f.first = read_indexer.get();
            f.second = 1;
        }
        else {
            excluded_frags.insert(*a);
        }

        /* This is not a multiread, so there must be at most one alignment. */
        assert(++a == AlignedReadIterator());
    }

    frag_counts.reserve(frag_counts.size() + frag_idx.size());
    for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
        frag_counts.push_back(f->second);
    }

    weight_matrix_entries.reserve(weight_matrix_entries.size() + num_entries);
    transcript_weights.reserve(transcript_weights.size() + locus->ts.size());
    for (TranscriptSetLocus::iterator t = locus->ts.begin();
            t != locus->ts.end(); ++t) {
        transcript_sequence_bias(*locus, *t);
        transcript_weights.push_back(std::make_pair(t->id, transcript_weight(*t)));

        for (MultireadSet::iterator r = multiread_set.begin();
             r != multiread_set.end(); ++r) {
            for (AlignedReadIterator a(*r->second); a != AlignedReadIterator(); ++a) {
                float w = fragment_weight(*t, *a);
                if (w > 0.0) {
                    multiread_index.push(r->first, Multiread::Cell(t->id, w, 1.0));
                }
            }
        }

        for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
            float w = fragment_weight(*t, f->first);
            if (w > 0.0) {
                weight_matrix_entries.push_back(
                        WeightMatrixEntry(t->id, f->second.first, w));
            }
        }
    }
}


void SamplerInitThread::transcript_sequence_bias(
                const SamplerInitInterval& locus,
                const Transcript& t)
{
    pos_t tlen = t.exonic_length();
    if ((size_t) tlen > mate1_seqbias[0].size()) {
        mate1_seqbias[0].resize(tlen);
        mate1_seqbias[1].resize(tlen);
        mate2_seqbias[0].resize(tlen);
        mate2_seqbias[1].resize(tlen);
    }

    if (fm.sb == NULL || locus.seq == NULL) {
        std::fill(mate1_seqbias[0].begin(), mate1_seqbias[0].begin() + tlen, 1.0);
        std::fill(mate1_seqbias[1].begin(), mate1_seqbias[1].begin() + tlen, 1.0);
        std::fill(mate2_seqbias[0].begin(), mate2_seqbias[0].begin() + tlen, 1.0);
        std::fill(mate2_seqbias[1].begin(), mate2_seqbias[1].begin() + tlen, 1.0);
        return;
    }

    t.get_sequence(tseq0, *locus.seq, fm.sb->getL(), fm.sb->getR());
    t.get_sequence(tseq1, *locus.seq, fm.sb->getR(), fm.sb->getL());
    tseq1.revcomp();

    for (pos_t pos = 0; pos < tlen; ++pos) {
        mate1_seqbias[0][pos] = fm.sb->get_mate1_bias(tseq0, pos + fm.sb->getL());
        mate1_seqbias[1][pos] = fm.sb->get_mate1_bias(tseq1, pos + fm.sb->getL());
        mate1_seqbias[0][pos] = fm.sb->get_mate2_bias(tseq0, pos + fm.sb->getL());
        mate1_seqbias[1][pos] = fm.sb->get_mate2_bias(tseq1, pos + fm.sb->getL());
    }
    std::reverse(mate1_seqbias[1].begin(), mate1_seqbias[1].begin() + tlen);
    std::reverse(mate2_seqbias[1].begin(), mate2_seqbias[1].begin() + tlen);
}



float SamplerInitThread::transcript_weight(const Transcript& t)
{
    pos_t trans_len = t.exonic_length();
    if ((size_t) trans_len + 1 > ws.size()) ws.resize(trans_len + 1);

    /* Set ws[k] to be the the number of fragmens of length k, weighted by
     * sequence bias. */
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = fm.frag_len_p(frag_len);

        /* Don't bother considering sequence bias if the fragment length
         * probability is extremely small (i.e., so that it suffocates any
         * effect from bias.). */
        if (frag_len_pr < constants::min_frag_len_pr) {
            ws[frag_len] = (float) (trans_len - frag_len + 1);
            continue;
        }

        /* TODO: The following logic assumes the library type is FR. We need
         * seperate cases to properly handle other library types, that
         * presumably exist somewhere. */

        ws[frag_len] = 0.0;

        float sp;
        sp = t.strand == strand_pos ? fm.strand_specificity :
                                      1.0 - fm.strand_specificity;
        for (pos_t pos = 0; pos <= trans_len - frag_len; ++pos) {
            ws[frag_len] += sp * mate1_seqbias[0][pos] * mate2_seqbias[1][pos + frag_len - 1];
        }

        sp = t.strand == strand_neg ? fm.strand_specificity :
                                      1.0 - fm.strand_specificity;
        for (pos_t pos = 0; pos <= trans_len - frag_len; ++pos) {
            ws[frag_len] += sp * mate2_seqbias[0][pos] * mate1_seqbias[1][pos + frag_len - 1];
        }
    }

    float w = 0.0;
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = fm.frag_len_p(frag_len);
        w += frag_len_pr * ws[frag_len];
    }

    return w;
}


float SamplerInitThread::fragment_weight(const Transcript& t,
                                         const AlignmentPair& a)
{
    pos_t frag_len = a.frag_len(t);
    pos_t trans_len = t.exonic_length();
    if (frag_len < 0.0) return 0.0;
    else if (frag_len == 0.0) {
        frag_len = std::min(trans_len, (pos_t) round(fm.frag_len_med()));
    }

    float w = fm.frag_len_p(frag_len);

    if (a.mate1) {
        pos_t offset = t.get_offset(a.mate1->strand == strand_pos ?
                                    a.mate1->start : a.mate1->end);
        assert(0 <= offset && offset < trans_len);
        w *= mate1_seqbias[a.mate1->strand][offset];
    }

    if (a.mate2) {
        pos_t offset = t.get_offset(a.mate2->strand == strand_pos ?
                                    a.mate2->start : a.mate2->end);
        assert(0 <= offset && offset < trans_len);
        w *= mate2_seqbias[a.mate2->strand][offset];
    }

    return w;
}


/* This class merges vectors of weight matrix entries generated by
 * initialization threads. */
class WeightMatrixEntryStream : public SparseMatEntryStream
{
    private:
        friend class HeapCmp;
        struct HeapCmp
        {
            HeapCmp(const WeightMatrixEntryStream& stream)
                : stream(stream)
            {}

            bool operator () (size_t i, size_t j)
            {
                if (stream.is[i] == stream.sources[i]->end()) return false;
                else if (stream.is[j] == stream.sources[j]->end()) return true;
                else {
                    if (stream.is[i]->i != stream.is[j]->i) {
                        return stream.is[i]->i < stream.is[j]->i;
                    }
                    else return stream.is[i]->j < stream.is[j]->j;
                }
            }

            const WeightMatrixEntryStream& stream;
        };

    public:
        WeightMatrixEntryStream(
                const std::vector<std::vector<WeightMatrixEntry>*>& sources)
            : cmp(*this)
            , is(sources.size())
            , sources(sources)
        {
            for (size_t i = 0; i < sources.size(); ++i) {
                is[i] = sources[i]->begin();
                if (is[i] != sources[i]->end()) {
                    h.push_back(i);
                }
            }

            std::make_heap(h.begin(), h.end(), cmp);
            next();
        }

        bool finished()
        {
            return h.empty();
        }

        void next()
        {
            size_t i = h.front();
            std::pop_heap(h.begin(), h.end(), cmp);
            h.pop_back();
            cur_row = is[i]->i;
            cur_col = is[i]->j;
            cur_val = is[i]->w;

            is[i]++;
            if (is[i] != sources[i]->end()) {
                h.push_back(i);
                std::push_heap(h.begin(), h.end());
            }
        }

        unsigned int row()
        {
            return cur_row;
        }

        unsigned int col()
        {
            return cur_col;
        }

        float value()
        {
            return cur_val;
        }

    private:
        HeapCmp cmp;
        std::vector<size_t> h;
        std::vector<std::vector<WeightMatrixEntry>::const_iterator> is;
        const std::vector<std::vector<WeightMatrixEntry>*>& sources;

        unsigned int cur_row, cur_col;
        float cur_val;
};



static unsigned int disjset_find(unsigned int* ds, unsigned int i)
{
    if (ds[i] == i) {
        return i;
    } else {
        return ds[i] = disjset_find(ds, ds[i]);
    }
}


static void disjset_union(unsigned int* ds, unsigned int i, unsigned int j)
{
    unsigned int a = disjset_find(ds, i);
    unsigned int b = disjset_find(ds, j);
    ds[b] = a;
}



Sampler::Sampler(const char* bam_fn, const char* fa_fn,
                 TranscriptSet& ts,
                 FragmentModel& fm)
    : ts(ts)
    , fm(fm)
{
    /* Producer/consumer queue of intervals containing indexed reads to be
     * processed. */
    Queue<SamplerInitInterval*> q(100);

    /* Assign matrix indices to reads. */
    Indexer read_indexer;

    /* Collect intervals */
    std::vector<SamplerInitInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(fm, read_indexer, q));
        threads.back()->start();
    }

    Logger::info("%lu loci", (unsigned long) intervals.size());

    sam_scan(intervals, bam_fn, fa_fn, "Estimating fragment weights");

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    /* Merge multiread indexes. */
    MultireadIndex multiread_index;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        for (MultireadIndex::iterator j = threads[i]->multiread_index.begin();
                j != threads[i]->multiread_index.end(); ++j) {
            multiread_index.push(j->first, j->second.cells);
        }

        /* free some memory */
        threads[i]->multiread_index.clear();
    }

    /* Add weight matrix entries for each multiread. */
    std::vector<WeightMatrixEntry> entries;
    for (MultireadIndex::iterator i = multiread_index.begin();
            i != multiread_index.end(); ++i) {
        std::deque<Multiread::Cell>& cells = i->second.cells;
        std::deque<Multiread::Cell>::iterator j;
        float sum_weight = 0.0;
        float sum_align_pr = 0.0;
        for (j = cells.begin(); j != cells.end(); ++j) {
            sum_weight += j->frag_weight * j->align_pr;
            sum_align_pr += j->align_pr;
        }

        if (sum_align_pr == 0.0 || sum_weight == 0.0) continue;

        for (j = cells.begin(); j != cells.end(); ++j) {
            j->align_pr /= sum_align_pr;
        }

        std::sort(cells.begin(), cells.end(), Multiread::CellTranOrd());

        unsigned int frag_idx = read_indexer.get();

        /* Sum alignments on the same transcripts. */
        float w = 0.0;
        std::deque<Multiread::Cell>::iterator j0;
        for (j0 = j = cells.begin(); j != cells.end(); ++j) {
            if (j->transcript_idx != j0->transcript_idx) {
                if (w > 0.0) {
                    entries.push_back(
                            WeightMatrixEntry(j0->transcript_idx, frag_idx, w));
                }
                j0 = j;
                w = 0.0;
            }

            w += j->align_pr * j->frag_weight;
        }
        if (w > 0.0) {
            entries.push_back(
                    WeightMatrixEntry(j0->transcript_idx, frag_idx, w));
        }
    }
    multiread_index.clear();

    /* Copy everything into the same vector. */
    for (size_t i = 0; i < constants::num_threads; ++i) {
        entries.reserve(entries.size() + threads[i]->weight_matrix_entries.size());
        for (std::vector<WeightMatrixEntry>::iterator j = threads[i]->weight_matrix_entries.begin();
                j != threads[i]->weight_matrix_entries.end(); ++j) {
            if (j->w > 0.0) {
                entries.push_back(*j);
            }
        }
    }

    /* Reassign read indices, now that zero entries are left out. */
    std::sort(entries.begin(), entries.end(), WeightMatrixEntry::ColOrd());
    unsigned int col = 0;
    unsigned int last_j = 0;
    for (std::vector<WeightMatrixEntry>::iterator i = entries.begin();
            i != entries.end(); ++i) {
        if (i->j != last_j && col > 0) ++col;
        i->j = col;
    }

    /* Find connected components. */
    size_t n = ts.size();
    size_t m = col + 1;

    /* ds is disjoint set data structure, where ds[i] holds a pointer to item
     * i's parent, and ds[i] = i, when the node is the root. */
    unsigned int* ds = new unsigned int [n + m];
    for (size_t i = 0; i < n + m; ++i) ds[i] = i;

    for (std::vector<WeightMatrixEntry>::iterator i = entries.begin();
            i != entries.end(); ++i) {
        if (i->w > 0.0) {
            disjset_union(ds, i->i, n + i->j);
        }
    }

    for (size_t i = 0; i < n + m; ++i) disjset_find(ds, i);

    /* 1. Relabel components
     * 2. Reassign column indices based on components. How?
     *    We need to sort the entries according to component.
     *    The easiest way may be to just add a component field to each entry.
     * */



    /* We can produce */

    //boost::unordered_set<unsigned int> components;
    //for (size_t i = 0; i < n + m; ++i) {
        //components.insert(disjset_find(ds, i));
    //}


#if 0
    /* Calculate the number of non-zero entries in the weight matrix. */
    std::vector<size_t> rowlens(ts.size(), 0);
    size_t ncol = 0;
    size_t nnz = multiread_entries.size();
    for (size_t i = 0; i < constants::num_threads; ++i) {
        std::vector<WeightMatrixEntry>::iterator j;
        for (j = threads[i]->weight_matrix_entries.begin();
                j != threads[i]->weight_matrix_entries.end(); ++j) {
            if (j->w > 0.0) {
                ++nnz;
                ++rowlens[j->i];
            }
        }
    }

    Logger::info("Weight matrix entries: %lu", (unsigned long) nnz);
#endif

    /* Reassign fragment indices, removing that can't have non-zero weight. */
    /* Can we do this with the stream trick?
     *
     *
     * Once we free the multiread indexes, I don't think memory is a big issue.
     * Let's just make one long entry list. Sort it. Reassign column indicies.
     *
     * Then, build a graph, find the connected components.
     *
     * This is going to involve another rasignment of column indices.
     */


#if 0
    std::vector<std::vector<WeightMatrixEntry>*> sources;
    sources.push_back(&multiread_entries);
    for (size_t i = 0; i < constants::num_threads; ++i) {
        sources.push_back(&threads[i]->weight_matrix_entries);
    }

    WeightMatrixEntryStream entry_stream(sources);
    W = new SparseMat(ts.size(), ncol, rowlens, entry_stream);
#endif


    /* Current problem:
     * What can be done about fragments with all zero weights.
     *
     * 1. Try to remove them here, which means reassigning indexs, etc.
     * 2. ...
     *
     * Fuck, we are going to have to do that anyway, since we removed a bunch of
     * reads when filtering the multireads most likely.
     */
}


Sampler::~Sampler()
{
}



