
/* Plan of attack:
 *
 * X 1. Make the indexer for reads start at multireads.size().
 *   2. Process MultiReadEntries into WeightEntries within each locus.
 *   3. In Sampler::Sampler figure out how to reassign indices to reads.
 *      We're likely going to have to do this twice.
 *      First, to remove zeros, and again to rearrange accorditing to component.
 */

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include "constants.hpp"
#include "hat-trie/hat-trie.h"
#include "linalg.hpp"
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



/* Safe c-style allocation. */

template <typename T>
T* malloc_or_die(size_t n)
{
    T* xs = reinterpret_cast<T*>(malloc(n * sizeof(T)));
    if (xs == NULL) {
        Logger::abort("Could not allocated %lu bytes.", n * sizeof(T));
    }
    return xs;
}


template <typename T>
T* realloc_or_die(T* xs, size_t n)
{
    T* ys = reinterpret_cast<T*>(realloc(xs, n * sizeof(T)));
    if (ys == NULL) {
        Logger::abort("Could not allocated %lu bytes.", n * sizeof(T));
    }
    return ys;
}



/* A row compressed sparse matrix.
 */
class WeightMatrix
{
    public:
        WeightMatrix(unsigned int nrow)
        {
            this->nrow = nrow;
            ncol = 0;

            rowlens = new unsigned int [nrow];
            std::fill(rowlens, rowlens+ nrow, 0);

            reserved = new unsigned int [nrow];
            std::fill(reserved, reserved + nrow, 0);

            rows = new float* [nrow];
            std::fill(rows, rows + nrow, (float*) NULL);

            idxs = new unsigned int* [nrow];
            std::fill(idxs, idxs + nrow, (unsigned int*) NULL);

            compacted = false;
        }

        ~WeightMatrix()
        {
            delete [] rowlens;
            delete [] reserved;
            if (compacted) {
                for (unsigned int i = 0; i < nrow; ++i) {
                    afree(rows[i]);
                    afree(idxs[i]);
                }
            }
            else {
                for (unsigned int i = 0; i < nrow; ++i) {
                    free(rows[i]);
                    free(idxs[i]);
                }
            }
            delete [] rows;
            delete [] idxs;
        }

        void push(unsigned int i, unsigned int j, float w)
        {
            if (rowlens[i] >= reserved[i]) {
                /* We use a pretty non-agressive resizing strategy, since there
                 * are potentially many arrays here and we want to avoid wasting
                 * much space.  */
                unsigned int newsize;
                if (reserved[i] == 0) {
                    newsize = 1;
                }
                else if (reserved[i] < 100) {
                    newsize = 2 * reserved[i];
                }
                else {
                    newsize = reserved[i] + 100;
                }

                rows[i] = realloc_or_die(rows[i], newsize);
                idxs[i] = realloc_or_die(idxs[i], newsize);
                reserved[i] = newsize;
            }

            idxs[i][rowlens[i]] = j;
            rows[i][rowlens[i]] = w;
            ++rowlens[i];
        };


        /* This function does two things: sorts the columns by index and
         * reallocates each array to be of exact size (to preserve space) and to
         * be 16-bytes aligned. It also reassigns column indices to remove empty
         * columns. */
        void compact()
        {
            /* Reallocate each row array */
            for (unsigned int i = 0; i < nrow; ++i) {
                unsigned int m = rowlens[i];
                if (m == 0) continue;

                float* newrow = reinterpret_cast<float*>(
                        aalloc(m * sizeof(float)));
                memcpy(newrow, rows[i], m * sizeof(float));
                free(rows[i]);
                rows[i] = newrow;

                unsigned int* newidx = reinterpret_cast<unsigned int*>(
                        aalloc(m * sizeof(unsigned int)));
                memcpy(newidx, idxs[i], m * sizeof(unsigned int));
                free(idxs[i]);
                idxs[i] = newidx;

                reserved[i] = m;
            }

            /* TODO: just testing sortrow here. */
            /* Sort each row by column */
            for (unsigned int i = 0; i < nrow; ++i) {
                sortrow(i);
            }

            /* Mark observed columns */
            ncol = 0;
            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    if (idxs[i][j] >= ncol) ncol = idxs[i][j] + 1;
                }
            }

            unsigned int* newidx = new unsigned int [ncol];
            memset(newidx, 0, ncol * sizeof(unsigned int));

            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    newidx[idxs[i][j]] = 1;
                }
            }

            /* Reassign column indexes */
            for (unsigned int i = 0, j = 0; i < ncol; ++i) {
                if (newidx[i]) newidx[i] = j++;
                else newidx[i] = j;
            }
            ncol = newidx[ncol - 1] + 1;

            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    idxs[i][j] = newidx[idxs[i][j]];
                }
            }

            delete [] newidx;

            /* Sort each row by column */
            for (unsigned int i = 0; i < nrow; ++i) {
                sortrow(i);
            }

            compacted = true;
        }

        // TODO: reorder_columns


        unsigned int nrow, ncol;

    private:
        unsigned int median_of_three(const unsigned int* idx,
                                     unsigned int a,
                                     unsigned int b,
                                     unsigned int c) const
        {
            unsigned int abc[3] = {a, b, c};
            if (idx[abc[0]] < idx[abc[1]]) std::swap(abc[0], abc[1]);
            if (idx[abc[1]] < idx[abc[2]]) std::swap(abc[1], abc[2]);
            if (idx[abc[0]] < idx[abc[1]]) std::swap(abc[0], abc[1]);
            return abc[1];
        }

        /* Sort a row by column index. A custom sort function is used since we
         * want to sort both rows[i] and idxs[i] by idxs[i] without any
         * intermediate steps. */
        void sortrow(unsigned int i)
        {
            unsigned int m = rowlens[i];
            if (m == 0) return;

            const unsigned int insert_sort_cutoff = 7;
            typedef std::pair<unsigned int, unsigned int> Range;
            typedef std::deque<Range> RangeStack;
            RangeStack s;
            float* row = rows[i];
            unsigned int* idx = idxs[i];

            s.push_back(std::make_pair(0, m));
            while (!s.empty()) {
                Range r = s.back();
                s.pop_back();

                unsigned int u = r.first;
                unsigned int v = r.second;

                if (u == v) continue;
                else if (v - u < insert_sort_cutoff) {
                    unsigned int a, minpos, minval;
                    while (u < v - 1) {
                        minpos = u;
                        minval = idx[u];
                        for (a = u + 1; a < v; ++a) {
                            if (idx[a] < minval) {
                                minpos = a;
                                minval = idx[a];
                            }
                        }
                        std::swap(idx[u], idx[minpos]);
                        std::swap(row[u], row[minpos]);
                        ++u;
                    }
                    continue;
                }

                unsigned int p = median_of_three(
                        idx, u, u + (v - u) / 2, v - 1);

                unsigned int pval = idx[p];

                std::swap(idx[p], idx[v - 1]);
                std::swap(row[p], row[v - 1]);

                unsigned int a, b;
                for (a = b = u; b < v - 1; ++b) {
                    if (idx[b] <= pval) {
                        std::swap(idx[a], idx[b]);
                        std::swap(row[a], row[b]);
                        ++a;
                    }
                }

                std::swap(idx[a], idx[v - 1]);
                std::swap(row[a], row[v - 1]);

                s.push_back(std::make_pair(u, a));
                s.push_back(std::make_pair(a + 1, v));
            }

            /* TODO: remove when we are sure this is working. */
            for (unsigned int j = 0; j < m - 1; ++j) {
                if (idx[j] >= idx[j + 1]) {
                    Logger::error("WeightMatrix::sortrow is broken.");
                }
            }
        }

        /* Number of entries in each row. */
        unsigned int* rowlens;

        /* Entries reserved in each row. */
        unsigned int* reserved;

        /* Row data. */
        float** rows;

        /* Columns indexed for each row. */
        unsigned int** idxs;

        /* Have arrays been resized to fit with aligned blocks of memory. */
        bool compacted;

        friend class WeightMatrixIterator;
};


struct WeightMatrixEntry
{
    unsigned int i, j;
    float w;
};


class WeightMatrixIterator :
    public boost::iterator_facade<WeightMatrixIterator,
                                  const WeightMatrixEntry,
                                  boost::forward_traversal_tag>
{
    public:
        WeightMatrixIterator()
            : owner(NULL)
        {

        }

        WeightMatrixIterator(const WeightMatrix& owner)
            : owner(&owner)
        {
            i = 0;
            k = 0;
        }


    private:
        friend class boost::iterator_core_access;

        void increment()
        {
            if (owner == NULL || i >= owner->nrow) return;

            ++k;
            while (i < owner->nrow && k >= owner->rowlens[i]) {
                ++i;
                k = 0;
            }

            entry.i = i;
            entry.j = owner->idxs[i][k];
            entry.w = owner->rows[i][k];
        }

        bool equal(const WeightMatrixIterator& other) const
        {
            if (owner == NULL || i >= owner->nrow) {
                return other.owner == NULL || other.i >= other.owner->nrow;
            }
            else if (other.owner == NULL || other.i >= other.owner->nrow) {
                return false;
            }
            else {
                return owner == other.owner &&
                       i == other.i &&
                       k == other.k;
            }
        }

        const WeightMatrixEntry& dereference() const
        {
            return entry;
        }

        WeightMatrixEntry entry;
        const WeightMatrix* owner;
        unsigned int i, k;
};


/* A trivial class to serve unique indices to multile threads. */
class Indexer
{
    public:
        Indexer(unsigned int first = 0)
            : next(first)
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
#if 0
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
#endif


struct MultireadEntry
{
    MultireadEntry(unsigned int multiread_num,
                   unsigned int transcript_idx,
                   float frag_weight,
                   float align_pr)
        : multiread_num(multiread_num)
        , transcript_idx(transcript_idx)
        , frag_weight(frag_weight)
        , align_pr(align_pr)
    {}


    bool operator < (const MultireadEntry& other) const
    {
        if (multiread_num != other.multiread_num) {
            return multiread_num < other.multiread_num;
        }
        else return transcript_idx < other.transcript_idx;
    }

    unsigned int multiread_num;
    unsigned int transcript_idx;
    float frag_weight;
    float align_pr;
};


/* Thread safe vector, allowing multiple threads to write.
 * In particular: modification through push_back and reserve_extra are safe. */
template <typename T>
class TSVec : public std::deque<T>
{
    public:
        void push_back(const T& x)
        {
            boost::lock_guard<boost::mutex> lock(mut);
            std::deque<T>::push_back(x);
        }

    private:
        boost::mutex mut;
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
                WeightMatrix& weight_matrix,
                Queue<SamplerInitInterval*>& q)
            : weight_matrix(weight_matrix)
            , fm(fm)
            , read_indexer(read_indexer)
            , q(q)
            , thread(NULL)
        {
        }

        ~SamplerInitThread()
        {
            /* Note: we are not free weight_matrix_entries and
             * multiread_entries. This get's done in Sampler::Sampler.
             * It's all part of the delicate dance involved it minimizing
             * memory use. */
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

        WeightMatrix& weight_matrix;

        typedef std::pair<unsigned int, unsigned int> FragIdxCount;
        std::vector<FragIdxCount> frag_counts;

        typedef std::pair<unsigned int, float> TransIdxCount;
        std::vector<TransIdxCount> transcript_weights;

        FragmentModel& fm;
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
                break;
            }
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

    std::vector<MultireadEntry> multiread_entries;
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
                    /* TODO: estimate alignment probability */
                    multiread_entries.push_back(
                            MultireadEntry(r->first, t->id, w, 1.0));
                }
            }
        }

        for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
            float w = fragment_weight(*t, f->first);
            if (w > 0.0) {
                weight_matrix.push(t->id, f->second.first, w);
            }
        }
    }

    /* Process multiread entries into weight matrix entries */
#if 0
    std::sort(multiread_entries.begin(), multiread_entries.end());
    for (std::vector<MultireadEntry>::iterator i = multiread_entries.begin();
            i != multiread_entries.end(); ++i) {
        std::vector<MultireadEntry>::iterator j, k = i + 1;
        while (k != multiread_entries.end() &&
                i->multiread_num == k->multiread_num) {
            ++k;
        }

        float sum_weight = 0.0;
        float sum_align_pr = 0.0;
        for (j = i; j != k; ++j) {
            sum_weight += j->frag_weight * j->align_pr;
            sum_align_pr += j->align_pr;
        }

        if (sum_align_pr == 0.0 || sum_weight == 0.0) continue;

        /* Sum weight of alignments on the same transcript. */
        float w = 0.0;
        std::vector<MultireadEntry>::iterator j0;
        for (j0 = j = i; j != k; ++j) {
            if (j->transcript_idx != j0->transcript_idx) {
                weight_matrix.push(j->transcript_idx, j0->multiread_num, w);
                j0 = j;
                w = 0.0;
            }

            w += j->align_pr * j->frag_weight;
        }
        if (w > 0.0) {
            weight_matrix.push(j0->transcript_idx, j0->multiread_num, w);
        }
    }
#endif
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
    Indexer read_indexer(fm.multireads.size());

    /* Collect intervals */
    std::vector<SamplerInitInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    WeightMatrix weight_matrix(ts.size());;

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(
                    fm,
                    read_indexer,
                    weight_matrix,
                    q));
        threads.back()->start();
    }

    Logger::info("Loci: %lu", (unsigned long) intervals.size());

    sam_scan(intervals, bam_fn, fa_fn, "Estimating fragment weights");

    /* Free a little space. */
    fm.multireads.clear();

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    Logger::info("Pre-compact columns: %lu",
            (unsigned long) weight_matrix.ncol);
    weight_matrix.compact();
    Logger::info("Post-compact columns: %lu",
            (unsigned long) weight_matrix.ncol);



    /* Find connected components. */

    unsigned int N = weight_matrix.ncol + weight_matrix.nrow;

    /* ds is disjoint set data structure, where ds[i] holds a pointer to item
     * i's parent, and ds[i] = i, when the node is the root. */
    unsigned int* ds = new unsigned int [N];
    for (size_t i = 0; i < N; ++i) ds[i] = i;

    /* TODO: we need to be able to iterate over matrix entries */
#if 0
    for (TSVec<WeightMatrixEntry>::iterator i = weight_matrix_entries.begin();
            i != weight_matrix_entries.end(); ++i) {
        if (i->w > 0.0) {
            disjset_union(ds, i->i, n + i->j);
        }
    }

    for (size_t i = 0; i < n + m; ++i) disjset_find(ds, i);

    /* Label components sequentially. */
    boost::unordered_map<unsigned int, unsigned int> component_label;
    for (size_t i = 0; i < n + m; ++i) {
        boost::unordered_map<unsigned int, unsigned int>::iterator c;
        c = component_label.find(ds[i]);
        if (c == component_label.end()) {
            component_label.insert(std::make_pair(ds[i], component_label.size()));
        }
    }
    Logger::info("Components: %lu", (unsigned long) component_label.size());
#endif



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
    //

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



