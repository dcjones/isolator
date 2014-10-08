
#include <boost/foreach.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <climits>

#include "constants.hpp"
#include "fastmath.hpp"
#include "hat-trie/hat-trie.h"
#include "logger.hpp"
#include "nlopt/nlopt.h"
#include "queue.hpp"
#include "read_set.hpp"
#include "sampler.hpp"
#include "samtools/sam.h"
#include "samtools/samtools_extra.h"
#include "shredder.hpp"

extern "C" {
#include "samtools/khash.h"
KHASH_MAP_INIT_STR(s, int)
static void supress_khash_unused_function_warnings() __attribute__ ((unused));
static void supress_khash_unused_function_warnings()
{
    UNUSED(kh_init_s);
    UNUSED(kh_destroy_s);
    UNUSED(kh_put_s);
    UNUSED(kh_clear_s);
    UNUSED(kh_del_s);
}
}


static double sq(double x)
{
    return x * x;
}


static double gamma_lnpdf(double alpha, double beta, double x)
{
    return alpha * fastlog(beta) -
           lgamma(alpha) +
           (alpha - 1) * fastlog(x) -
           beta * x;
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


void assert_finite(double x)
{
    if (!boost::math::isfinite(x)) {
        Logger::abort("Unexpected non-finite value.");
    }
}


struct MultireadFrag
{
    MultireadFrag(unsigned int multiread_num,
                  unsigned int frag_idx,
                  float align_pr)
        : multiread_num(multiread_num)
        , frag_idx(frag_idx)
        , align_pr(align_pr)
    {
    }

    unsigned int multiread_num;
    unsigned int frag_idx;
    float align_pr;

    bool operator < (const MultireadFrag& other) const
    {
        if (multiread_num != other.multiread_num) {
            return multiread_num < other.multiread_num;
        }
        else return frag_idx < other.frag_idx;
    }
};


/* A (fragment indx, count) pair. */
typedef std::pair<unsigned int, unsigned int> FragIdxCount;


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
         * columns.
         *
         * Returns:
         *   A map of the previous column indexes onto new indexs. The caller is
         *   responsible for freeing this with delete [].
         *
         * */
        unsigned int* compact(size_t* oldncol)
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

            /* Mark observed columns */
            ncol = 0;
            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    if (idxs[i][j] >= ncol) ncol = idxs[i][j] + 1;
                }
            }
            if (ncol == 0) {
                Logger::abort("Dataset has no viable fragments. "
                              "This is most likely a bug.");
            }

            unsigned int* newidx = new unsigned int [ncol];
            memset(newidx, 0, ncol * sizeof(unsigned int));
            *oldncol = ncol;

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

            /* Sort each row by column */
            for (unsigned int i = 0; i < nrow; ++i) {
                sortrow(i);
            }

            compacted = true;

            return newidx;
        }


        /* Reorder columns given an array of the form map[i] = j, mapping column
         * i to column j. */
        void reorder_columns(const unsigned int* idxmap)
        {
            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    idxs[i][j] = idxmap[idxs[i][j]];
                }
                sortrow(i);
            }
        }

        unsigned int nrow, ncol;

        /* Number of entries in each row. */
        unsigned int* rowlens;

        /* Row data. */
        float** rows;

        /* Columns indexed for each row. */
        unsigned int** idxs;

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
        }

        /* Entries reserved in each row. */
        unsigned int* reserved;

        /* Have arrays been resized to fit with aligned blocks of memory. */
        bool compacted;

        friend class WeightMatrixIterator;
};


struct WeightMatrixEntry
{
    WeightMatrixEntry()
        : i(-1)
        , j(-1)
        , w(NAN)
    {
    }

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
            while (i < owner.nrow && 0 >= owner.rowlens[i]) {
                ++i;
            }

            if (i < owner.nrow) {
                entry.i = i;
                entry.j = owner.idxs[i][0];
                entry.w = owner.rows[i][0];
            }
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

            if (i < owner->nrow) {
                entry.i = i;
                entry.j = owner->idxs[i][k];
                entry.w = owner->rows[i][k];
            }
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


struct MultireadEntry
{
    MultireadEntry(unsigned int multiread_num,
                   unsigned int alignment_num,
                   unsigned int transcript_idx,
                   float frag_weight,
                   float align_pr)
        : multiread_num(multiread_num)
        , alignment_num(alignment_num)
        , transcript_idx(transcript_idx)
        , frag_weight(frag_weight)
        , align_pr(align_pr)
    {}


    bool operator < (const MultireadEntry& other) const
    {
        if (multiread_num != other.multiread_num) {
            return multiread_num < other.multiread_num;
        }
        else if (alignment_num != other.alignment_num) {
            return alignment_num < other.alignment_num;
        }
        else return transcript_idx < other.transcript_idx;
    }

    unsigned int multiread_num;
    unsigned int alignment_num;
    unsigned int transcript_idx;
    float frag_weight;
    float align_pr;
};


/* Thread safe vector, allowing multiple threads to write.
 * In particular: modification through push_back and reserve_extra are safe. */
template <typename T>
class TSVec : public std::vector<T>
{
    public:
        void push_back(const T& x)
        {
            boost::lock_guard<boost::mutex> lock(mut);
            std::vector<T>::push_back(x);
        }


        void reserve_extra(size_t n)
        {
            boost::lock_guard<boost::mutex> lock(mut);
            std::vector<T>::reserve(std::vector<T>::size() + n);
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
            start = std::max<pos_t>(0, ts.min_start - 1000);
            end = ts.max_end + 1000;
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
            else if (start != other.start) {
                return start < other.start;
            }
            else {
                return end < other.end;
            }
        }

        TranscriptSetLocus ts;
        ReadSet rs;
        boost::shared_ptr<twobitseq> seq;

        int32_t tid;
        pos_t start, end;

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
    while (samread(bam_f, b) >= 0) {
        ++read_num;

        if (read_num % 1000 == 0) {
            file_pos = samtell(bam_f);
            if (file_pos >= last_file_pos + input_block_size && input_size > 0) {
                Logger::get_task(task_name).inc();
                last_file_pos = file_pos;
            }
        }

        if (b->core.flag & BAM_FUNMAP || b->core.tid < 0) continue;
        if (b->core.mtid != -1 && b->core.tid != b->core.mtid) continue;

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

            for (j = j0; j < n &&intervals[j]->tid <= b->core.tid; ++j) {
                if (intervals[j]->tid < b->core.tid) continue;
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

            if (b->core.pos < intervals[j]->start) break;
            if (b->core.pos > intervals[j]->end) {
                if (j == j0) {
                    intervals[j0++]->finish();
                }
                continue;
            }

            pos_t b_end = (pos_t) bam_calend(&b->core, bam1_cigar(b)) - 1;
            if (b_end <= intervals[j]->end) {
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


class FragWeightEstimationThread
{
    public:
        FragWeightEstimationThread(
                FragmentModel& fm,
                Indexer& read_indexer,
                WeightMatrix& weight_matrix,
                TSVec<FragIdxCount>& frag_counts,
                TSVec<MultireadFrag>& multiread_frags,
                double* transcript_weights,
                double* transcript_gc,
                Queue<SamplerInitInterval*>& q)
            : weight_matrix(weight_matrix)
            , frag_counts(frag_counts)
            , multiread_frags(multiread_frags)
            , transcript_weights(transcript_weights)
            , transcript_gc(transcript_gc)
            , fm(fm)
            , read_indexer(read_indexer)
            , q(q)
            , thread(NULL)
            , seqbias_size(0)
            , posbias_size(0)
        {
            seqbias[0] = seqbias[1] = NULL;
            posbias = NULL;
            if (fm.frag_len_dist) {
                frag_len_dist = new EmpDist(*fm.frag_len_dist);
            }
            else {
                frag_len_dist = NULL;
            }
        }

        ~FragWeightEstimationThread()
        {
            delete frag_len_dist;
            afree(seqbias[0]);
            afree(seqbias[1]);
            afree(posbias);
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
            thread = new boost::thread(boost::bind(&FragWeightEstimationThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        double frag_len_p(pos_t frag_len)
        {
            if (frag_len_dist) return frag_len_dist->pdf(frag_len);
            else return fm.frag_len_p(frag_len);
        }

        double frag_len_c(pos_t frag_len)
        {
            if (frag_len_dist) return frag_len_dist->cdf(frag_len);
            else return fm.frag_len_c(frag_len);
        }

        void process_locus(SamplerInitInterval* locus);

        /* Compute sequence bias for both mates on both strand.
         *
         * Results are stored in seqbias. */
        void transcript_sequence_bias(const SamplerInitInterval& locus,
                                      const Transcript& t);

        /* Compute GC bias for a fixed fragment length.
         *
         * Results are stored in gcbias. */
        void transcript_gc_bias(
            const SamplerInitInterval& locus,
            const Transcript& t, pos_t frag_len);

        void transcript_tp_bias(
            const Transcript& t, pos_t frag_len);

        /* Compute the sum of the weights of all fragmens in the transcript.
         *
         * This assumes that transcript_sequence_bias has already been called,
         * and the seqbias arrays set for t. */
        float transcript_weight(const SamplerInitInterval& locus,
                                const Transcript& t);

        /* Compute (a number proportional to) the probability of observing the
         * given fragment from the given transcript, assuming every transcript
         * is equally expressed. */
        float fragment_weight(
                const Transcript& t,
                const AlignmentPair& a);

        WeightMatrix& weight_matrix;
        TSVec<FragIdxCount>& frag_counts;
        TSVec<MultireadFrag>& multiread_frags;
        double* transcript_weights;
        double* transcript_gc;

        FragmentModel& fm;
        Indexer& read_indexer;
        Queue<SamplerInitInterval*>& q;
        boost::thread* thread;

        /* Copy the fragment length distribution to avoid contention between
         * threads. */
        EmpDist* frag_len_dist;

        /* Temprorary space for computing sequence bias, indexed by
           sense (0) / antisense (1) */
        float* seqbias[2];
        size_t seqbias_size;

        // fragment gc bias and 3' bias
        float* posbias;
        size_t posbias_size;

        /* Exonic length of the transcript whos bias is stored in seqbias. */
        pos_t tlen;
        float tw;

        /* Temporary space for holding transcript sequences. */
        twobitseq tseq0; /* + strand */
        twobitseq tseq1; /* reverse complement of tseq0 */

        /* Temporary space used when computing transcript weight. */
        std::vector<float> ws;

        friend class Sampler;
};


void FragWeightEstimationThread::process_locus(SamplerInitInterval* locus)
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
        if (fm.blacklist.get(r->first) >= 0) continue;

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
            /* we assign a read index later so we can have indexes roughly
             * ordered by position. */
            f.first = 0;
            f.second = 1;
        }
        else {
            excluded_frags.insert(*a);
        }

        /* This is not a multiread, so there must be at most one alignment. */
        assert(++a == AlignedReadIterator());
    }

    /* Asign indexes */
    for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
        f->second.first = read_indexer.get();
    }

    for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
        /* fragment with count == 1 are implicit */
        if (f->second.second > 1) {
            frag_counts.push_back(f->second);
        }
    }

    std::vector<MultireadEntry> multiread_entries;
    for (TranscriptSetLocus::iterator t = locus->ts.begin();
            t != locus->ts.end(); ++t) {

        pos_t L = constants::seqbias_left_pos, R = constants::seqbias_right_pos;
        if (locus->seq) {
            t->get_sequence(tseq0, *locus->seq, L, R);
            t->get_sequence(tseq1, *locus->seq, R, L);
            tseq1.revcomp();
        }

        transcript_sequence_bias(*locus, *t);
        transcript_weights[t->id] = transcript_weight(*locus, *t);

        /* If a transcript is very unlikely to have been sequenced at all,
         * ignore it. Otherwise quite a few outliers will be generated in which
         * the annotated transcript is much shorter that what is actually being
         * transcribed. */
        if (transcript_weights[t->id] < constants::min_transcript_weight) {
            continue;
        }

        for (MultireadSet::iterator r = multiread_set.begin();
             r != multiread_set.end(); ++r) {
            unsigned int alignment_num = 0;
            for (AlignedReadIterator a(*r->second); a != AlignedReadIterator();
                    ++a, ++alignment_num) {
                float w = fragment_weight(*t, *a);
                if (w > 0.0) {
                    double align_pr = 1.0;

                    // TODO: we really need to do a better job of this shit
                    if (a->mate1) {
                        align_pr *= pow(0.01, a->mate1->mismatch_count);
                        //align_pr *= 1.0 - pow(10.0, -0.1 * a->mate1->mapq);
                    }

                    if (a->mate2) {
                        align_pr *= pow(0.01, a->mate2->mismatch_count);
                        //align_pr *= 1.0 - pow(10.0, -0.1 * a->mate2->mapq);
                    }

                    // avoid numerical issues
                    align_pr = std::max<double>(align_pr, 1e-6);

                    multiread_entries.push_back(
                            MultireadEntry(r->first, alignment_num, t->id, w, align_pr));
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
    std::sort(multiread_entries.begin(), multiread_entries.end());
    for (std::vector<MultireadEntry>::iterator k, j, i = multiread_entries.begin();
            i != multiread_entries.end(); i = k) {
        k = i + 1;
        while (k != multiread_entries.end() &&
               i->multiread_num == k->multiread_num &&
               i->alignment_num == k->alignment_num) {
            ++k;
        }

        float sum_weight = 0.0;
        float sum_align_pr = 0.0;
        for (j = i; j != k; ++j) {
            sum_weight += j->frag_weight;
            sum_align_pr += j->align_pr;
        }
        if (sum_weight == 0.0 || sum_align_pr == 0.0) {
            continue;
        }

        unsigned int frag_idx = read_indexer.get();
        multiread_frags.push_back(MultireadFrag(
                    i->multiread_num, frag_idx, sum_align_pr));

        /* Sum weight of alignments on the same transcript. */
        float w = 0.0;
        std::vector<MultireadEntry>::iterator j0;

        double align_pr_sum = 0.0;
        for (j0 = j = i; j != k; ++j) {
            if (j->transcript_idx != j0->transcript_idx) {
                if (w > 0.0) {
                    weight_matrix.push(j0->transcript_idx, frag_idx,
                            w / align_pr_sum);
                }
                j0 = j;
                w = 0.0;
                align_pr_sum = 0.0;
            }

            if (j->frag_weight > 0.0) {
                w += j->frag_weight;
                align_pr_sum += 1.0;
            }
        }
        if (w > 0.0) {
            weight_matrix.push(j0->transcript_idx, frag_idx,
                    w / align_pr_sum);
        }
    }
}


void FragWeightEstimationThread::transcript_sequence_bias(
                const SamplerInitInterval& locus,
                const Transcript& t)
{
    tlen = t.exonic_length();
    if ((size_t) tlen > seqbias_size) {
        afree(seqbias[0]);
        afree(seqbias[1]);
        seqbias[0] = reinterpret_cast<float*>(aalloc(tlen * sizeof(float)));
        seqbias[1] = reinterpret_cast<float*>(aalloc(tlen * sizeof(float)));
        seqbias_size = tlen;
    }

    std::fill(seqbias[0], seqbias[0] + tlen, 1.0);
    std::fill(seqbias[1], seqbias[1] + tlen, 1.0);

    if (fm.sb[1] == NULL || locus.seq == NULL){
        transcript_gc[t.id] = 0.5;
        return;
    }

    pos_t L = constants::seqbias_left_pos;
    transcript_gc[t.id] = tseq0.gc_count(L, L + tlen - 1) / (double) tlen;

    if (t.strand == strand_pos) {
        for (pos_t pos = 0; pos < tlen; ++pos) {
            seqbias[0][pos] = fm.sb[0]->get_bias(tseq0, pos + L);
            seqbias[1][pos] = fm.sb[1]->get_bias(tseq1, pos + L);
        }
        std::reverse(seqbias[1], seqbias[1] + tlen);
    }
    else {
        for (pos_t pos = 0; pos < tlen; ++pos) {
            seqbias[0][pos] = fm.sb[0]->get_bias(tseq1, pos + L);
            seqbias[1][pos] = fm.sb[1]->get_bias(tseq0, pos + L);
        }
        std::reverse(seqbias[0], seqbias[0] + tlen);
    }
}


void FragWeightEstimationThread::transcript_gc_bias(
    const SamplerInitInterval& locus, const Transcript& t, pos_t frag_len)
{
    if (!fm.gcbias || !locus.seq) {
        return;
    }

    tlen = t.exonic_length();
    pos_t L = constants::seqbias_left_pos;
    unsigned int gc_count = 0;
    for (pos_t pos = 0; pos < tlen; ++pos) {
        if (tseq0.isgc(pos + L)) ++gc_count;
        if (pos >= frag_len) {
            if (tseq0.isgc(pos - frag_len + L)) --gc_count;
        }
        if (pos >= frag_len - 1) {
            double gc = (double) gc_count / (double) frag_len;
            posbias[pos - frag_len + 1] *= fm.gcbias->get_bias(gc);
        }
    }
}


void FragWeightEstimationThread::transcript_tp_bias(const Transcript& t, pos_t frag_len)
{
    if (!fm.tpbias || fm.tpbias->p == 0.0) return;

    tlen = t.exonic_length();
    double p = fm.tpbias->p;
    double bias = fm.tpbias->get_bias(frag_len - 1);
    if (t.strand == strand_pos) {
        for (pos_t pos = tlen - frag_len; pos >= 0; --pos) {
            posbias[pos] *= bias;
            bias *= 1 - p;
        }
    }
    else {
        for (pos_t pos = frag_len - 1; pos < tlen; ++pos) {
            posbias[pos - frag_len + 1] *= bias;
            bias *= 1 - p;
        }
    }
}


float FragWeightEstimationThread::transcript_weight(
    const SamplerInitInterval& locus, const Transcript& t)
{
    pos_t trans_len = t.exonic_length();
    if ((size_t) trans_len + 1 > ws.size()) ws.resize(trans_len + 1);

    if ((size_t) tlen > posbias_size) {
        afree(posbias);
        posbias = reinterpret_cast<float*>(aalloc(tlen * sizeof(float)));
        posbias_size = tlen;
    }

    pos_t tlen = t.exonic_length();

    /* Set ws[k] to be the the number of fragmens of length k, weighted by
     * sequence bias. */
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = frag_len_p(frag_len);

        /* Don't bother considering sequence bias if the fragment length
         * probability is extremely small (i.e., so that it suffocates any
         * effect from bias.). */
        if (frag_len_pr < constants::transcript_len_min_frag_pr) {
            ws[frag_len] = (float) (trans_len - frag_len + 1);
            continue;
        }

        std::fill(posbias, posbias + tlen, 1.0);
        transcript_gc_bias(locus, t, frag_len);
        transcript_tp_bias(t, frag_len);

        /* TODO: The following logic assumes the library type is FR. We need
         * seperate cases to properly handle other library types, that
         * presumably exist somewhere. */

        ws[frag_len] = 0.0;

        ws[frag_len] +=
            fm.strand_specificity *
            dot(&seqbias[0][0],
                &seqbias[1][frag_len - 1],
                posbias,
                trans_len - frag_len + 1);

        ws[frag_len] +=
            (1.0 - fm.strand_specificity) *
            dot(&seqbias[1][0],
                &seqbias[0][frag_len - 1],
                posbias,
                trans_len - frag_len + 1);

        if (ws[frag_len] == 0.0) {
            Logger::abort("Numerical error in transcript weight computation.");
        }
    }

    tw = 0.0;
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = frag_len_p(frag_len) / frag_len_c(tlen);
        tw += frag_len_pr * ws[frag_len];
    }

    if (!boost::math::isfinite(tw) || tw <= constants::min_transcript_weight) {
        tw = constants::min_transcript_weight;
    }

    return tw;
}


float FragWeightEstimationThread::fragment_weight(const Transcript& t,
                                                  const AlignmentPair& a)
{
    pos_t frag_len = a.frag_len(t);
    pos_t L = constants::seqbias_left_pos;

    if (frag_len < 0) {
        return 0.0;
    }
    else if (frag_len == 0) {
        pos_t max_frag_len;
        const Alignment* mate = a.mate1 ? a.mate1 : a.mate2;

        if (mate->strand == strand_pos) {
            max_frag_len = t.max_end - mate->start + 1;
        }
        else {
            max_frag_len = mate->end - t.min_start + 1;
        }

        frag_len = std::min(max_frag_len, (pos_t) round(fm.frag_len_med()));
    }

    if (frag_len > tlen) {
        return 0.0;
    }

    float w = 1.0;
    if (a.mate1 && a.mate2) {
        pos_t offset1 = t.get_offset(a.mate1->strand == strand_pos ?
                                     a.mate1->start : a.mate1->end);
        if (offset1 < 0 || offset1 >= tlen) return 0.0;
        w *= seqbias[a.mate1->strand == t.strand ? 0 : 1][offset1];

        pos_t offset2 = t.get_offset(a.mate2->strand == strand_pos ?
                                     a.mate2->start : a.mate2->end);
        if (offset2 < 0 || offset2 >= tlen) return 0.0;
        w *= seqbias[a.mate2->strand == t.strand ? 0 : 1][offset2];

        pos_t offset = t.get_offset(
                std::min<pos_t>(a.mate1->start, a.mate2->start));
        if (fm.gcbias) {
            float gc = tseq0.gc_count(
                    L + offset, L + offset + frag_len - 1) / (float) frag_len;
            w *= fm.gcbias->get_bias(gc / frag_len);
        }

    }
    else {
        if (a.mate1) {
            pos_t offset = t.get_offset(a.mate1->strand == strand_pos ?
                                        a.mate1->start : a.mate1->end);
            if (offset < 0 || offset >= tlen) return 0.0;
            w *= seqbias[a.mate1->strand == t.strand ? 0 : 1][offset];

            if (fm.gcbias) {
                if (a.mate1->strand == strand_neg) {
                    offset = std::max<pos_t>(0, offset - frag_len + 1);
                }
                float gc = tseq0.gc_count(
                        L + offset,
                        std::min<pos_t>(tlen - 1, L + offset + frag_len - 1));
                w *= fm.gcbias->get_bias(gc / frag_len);
            }
        }

        if (a.mate2) {
            pos_t offset = t.get_offset(a.mate2->strand == strand_pos ?
                                        a.mate2->start : a.mate2->end);
            if (offset < 0 || offset >= tlen) return 0.0;
            w *= seqbias[a.mate2->strand == t.strand ? 0 : 1][offset];
            assert_finite(w);

            if (fm.gcbias) {
                if (a.mate2->strand == strand_neg) {
                    offset = std::max<pos_t>(0, offset - frag_len + 1);
                }
                float gc = tseq0.gc_count(
                        L + offset,
                        std::min<pos_t>(tlen - 1, L + offset + frag_len - 1));
                w *= fm.gcbias->get_bias(gc / frag_len);
                assert_finite(w);
            }
        }
    }

    // 3' bias
    if (t.strand == strand_pos && fm.tpbias && fm.tpbias->p != 0.0) {
        pos_t tpos;
        if (a.mate1 && a.mate1->strand == 0) {
            tpos = t.get_offset(a.mate1->start);
        }
        else if (a.mate2 && a.mate2->strand == 0) {
            tpos = t.get_offset(a.mate2->start);
        }
        else {
            const Alignment* mate = a.mate1 ? a.mate1 : a.mate2;
            tpos = std::max<pos_t>(0, t.get_offset(mate->end) - frag_len);
        }
        assert(0 <= tpos && tpos < tlen);
        w *= fm.tpbias->get_bias(tlen - tpos - 1);
        assert_finite(w);
    }
    else if (t.strand == strand_neg && fm.tpbias && fm.tpbias->p != 0.0) {
        pos_t tpos;
        if (a.mate1 && a.mate1->strand == 1) {
            tpos = t.get_offset(a.mate1->end);
        }
        else if (a.mate2 && a.mate2->strand == 1) {
            tpos = t.get_offset(a.mate2->end);
        }
        else {
            const Alignment* mate = a.mate1 ? a.mate1 : a.mate2;
            tpos = std::min<pos_t>(tlen - 1, t.get_offset(mate->start) + frag_len);
        }
        assert(0 <= tpos && tpos < tlen);
        w *= fm.tpbias->get_bias(tpos);
        assert_finite(w);
    }

    // strand-specificity
    if (a.mate1) {
        w *= a.mate1->strand == t.strand ?
            fm.strand_specificity : (1.0 - fm.strand_specificity);
    }
    else {
        w *= a.mate2->strand != t.strand ?
            fm.strand_specificity : (1.0 - fm.strand_specificity);
    }

    float frag_len_pr = frag_len_p(frag_len);
    if (frag_len_pr < constants::min_frag_len_pr) {
        return 0.0;
    }
    frag_len_pr /= frag_len_c(tlen);

    if (ws[frag_len] == 0.0 || frag_len_pr * w < constants::min_frag_weight) {
        return 0.0;
    }

    w = frag_len_pr * w / ws[frag_len];

    if (!boost::math::isfinite(w)) {
        Logger::abort("Non-finite fragment weight computed. %e",
                      frag_len_pr);
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


/* A comparison operator used in Sampler::Sampler. */
struct ComponentCmp
{
    ComponentCmp(unsigned int* ds)
        : ds(ds)
    {
    }

    bool operator () (unsigned int i, unsigned int j) {
        return ds[i] < ds[j];
    }

    unsigned int* ds;
};


class InterTgroupSampler : public Shredder
{
    public:
        InterTgroupSampler(Sampler& S)
            : Shredder(1e-6, 1.0 - 1e-6, 1e-6)
            , S(S)
        {

        }

        ~InterTgroupSampler()
        {
        }

        double sample(rng_t& rng, unsigned int c, unsigned int u, unsigned int v,
                      bool optimize_state)
        {
            double x0 = S.tgroupmix[u] / (S.tgroupmix[u] + S.tgroupmix[v]);

            if (upper_limit <= lower_limit) return x0;

            unsigned int component_size = S.component_frag[c + 1] - S.component_frag[c];

            // find the range of fragments affected by changes in x
            unsigned int f0 = S.component_frag[c+1];
            unsigned int f1 = S.component_frag[c];
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                if (S.weight_matrix->rowlens[tid] > 0) {
                    f0 = std::min(f0, S.weight_matrix->idxs[tid][0]);
                    f1 = std::max(f1, 1 + S.weight_matrix->idxs[tid][S.weight_matrix->rowlens[tid] - 1]);
                }
            }

            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                if (S.weight_matrix->rowlens[tid] > 0) {
                    f0 = std::min(f0, S.weight_matrix->idxs[tid][0]);
                    f1 = std::max(f1, 1 + S.weight_matrix->idxs[tid][S.weight_matrix->rowlens[tid] - 1]);
                }
            }
            if (f0 > f1) f0 = f1 = S.component_frag[c];
            f0 -= S.component_frag[c];
            f1 -= S.component_frag[c];

            // Round f0 down to the nearst 32-byte boundry so we can use avx/sse.
            f0 &= 0xfffffff8;

            this->x0 = x0;
            this->lpf01 = dotlog(S.frag_counts[c] + f0, S.frag_probs[c] + f0, f1 - f0) / M_LOG2E;
            this->c = c;
            this->u = u;
            this->v = v;
            this->f0 = f0;
            this->f1 = f1;
            this->tgroupmix_uv = S.tgroupmix[u] + S.tgroupmix[v];
            this->lp0 = dotlog(S.frag_counts[c], S.frag_probs[c], component_size) / M_LOG2E;

            prior_lp0 = 0.0;
            if (S.use_priors) {
                double tgroup_xu = fastlog(S.tgroupmix[u] * S.cmix[c]) -
                                   fastlog(S.tgroup_scaling[u]);
                prior_lp0 += tgroup_prior.f(S.hp.tgroup_mu[u], S.hp.tgroup_sigma[u],
                                            &tgroup_xu, 1);
                double tgroup_xv = fastlog(S.tgroupmix[v] * S.cmix[c]) -
                                   fastlog(S.tgroup_scaling[v]);
                prior_lp0 += tgroup_prior.f(S.hp.tgroup_mu[v], S.hp.tgroup_sigma[v],
                                            &tgroup_xv, 1);
                this->lp0 += prior_lp0;

                // jacobian
                prior_lp0 -= fastlog(x0);
                prior_lp0 -= fastlog(1 - x0);
            }

            if (optimize_state) {
                return Shredder::optimize(x0);
            }
            else {
                return Shredder::sample(rng, x0);
            }
        }

    protected:
        double f(double x, double& d)
        {
            double tgroupmix_u = x * tgroupmix_uv;
            double tgroupmix_v = (1 - x) * tgroupmix_uv;

            double u_delta = tgroupmix_u / S.tgroupmix[u];
            double v_delta = tgroupmix_v / S.tgroupmix[v];

            acopy(S.frag_probs_prop[c] + f0, S.frag_probs[c] + f0, (f1 - f0) * sizeof(float));
            d = 0.0;

            // recompute fragment probabilities
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                asxpy(S.frag_probs_prop[c],
                      S.weight_matrix->rows[tid],
                      S.tmix[tid] * u_delta - S.tmix[tid],
                      S.weight_matrix->idxs[tid],
                      S.component_frag[c],
                      S.weight_matrix->rowlens[tid]);
            }

            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                asxpy(S.frag_probs_prop[c],
                      S.weight_matrix->rows[tid],
                      S.tmix[tid] * v_delta - S.tmix[tid],
                      S.weight_matrix->idxs[tid],
                      S.component_frag[c],
                      S.weight_matrix->rowlens[tid]);
            }

            // gradient
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                const unsigned int rowlen = S.weight_matrix->rowlens[tid];
                const double z =
                    S.tgroup_tmix[tid] *
                    tgroupmix_uv;
                const unsigned int frag_offset = S.component_frag[c];
                const unsigned int* idxs = S.weight_matrix->idxs[tid];

                d += z * asxtydsz(S.frag_counts[c], S.weight_matrix->rows[tid],
                                  S.frag_probs_prop[c], idxs, frag_offset, rowlen);
            }

            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                const unsigned int rowlen = S.weight_matrix->rowlens[tid];
                const double z =
                    S.tgroup_tmix[tid] *
                    tgroupmix_uv;
                const unsigned int frag_offset = S.component_frag[c];
                const unsigned int* idxs = S.weight_matrix->idxs[tid];

                d -= z * asxtydsz(S.frag_counts[c], S.weight_matrix->rows[tid],
                                  S.frag_probs_prop[c], idxs, frag_offset, rowlen);
            }

            // prior probability
            double prior_lp_delta = 0.0;
            if (S.use_priors) {
                prior_lp_delta -= prior_lp0;
                double xu = S.cmix[c] * tgroupmix_u;
                double logxu = fastlog(xu) - fastlog(S.tgroup_scaling[u]);
                double mu_u = S.hp.tgroup_mu[u];
                prior_lp_delta += tgroup_prior.f(mu_u, S.hp.tgroup_sigma[u], &logxu, 1);
                prior_lp_delta -= fastlog(x); // jacobian

                double xv = S.cmix[c] * tgroupmix_v;
                double logxv = fastlog(xv) - fastlog(S.tgroup_scaling[v]);
                double mu_v = S.hp.tgroup_mu[v];
                prior_lp_delta += tgroup_prior.f(mu_v, S.hp.tgroup_sigma[v], &logxv, 1);
                prior_lp_delta -= fastlog(1 - x); // jacobian

                // derivative of normal log-pdf
                d += (mu_u - logxu) / (sq(S.hp.tgroup_sigma[u]) * x);
                d += (mu_v - logxv) / (sq(S.hp.tgroup_sigma[v]) * (x - 1));

                // derivative of jacobian
                d -= 1/x;
                d -= 1/(x-1);
            }

            return lp0 - lpf01 +
                   dotlog(S.frag_counts[c] + f0, S.frag_probs_prop[c] + f0, f1 - f0) / M_LOG2E +
                   prior_lp_delta;
        }

    private:
        NormalLogPdf tgroup_prior;

        double x0;

        // log-probability evaluated at x0
        double lp0;

        // log-probability of the fragments in [f0, f1]
        double lpf01;
        unsigned int c;
        unsigned int u;
        unsigned int v;
        unsigned int f0, f1;
        double tgroupmix_uv;
        double prior_lp0;

        Sampler& S;
};


double intertranscriptsampler_nlopt_objective(unsigned int _n, const double* _x,
                                              double* _grad, void* data);


class InterTranscriptSampler
{
    public:
        InterTranscriptSampler(Sampler& S)
            : opt(NULL)
            , S(S)
            , lower_limit(constants::zero_eps)
            , upper_limit(1.0 - constants::zero_eps)
        {
            //opt = nlopt_create(NLOPT_LD_SLSQP, 1);
            opt = nlopt_create(NLOPT_LN_SBPLX, 1);
            nlopt_set_lower_bounds(opt, &lower_limit);
            nlopt_set_upper_bounds(opt, &upper_limit);
            nlopt_set_max_objective(opt, intertranscriptsampler_nlopt_objective,
                                    reinterpret_cast<void*>(this));

            nlopt_set_ftol_abs(opt, 1e-2);
        }

        ~InterTranscriptSampler()
        {
            nlopt_destroy(opt);
        }

        double sample(unsigned int u, unsigned int v)
        {
            prepare_sample_optimize(u, v);
            double x0 = S.tmix[u] / (S.tmix[u] + S.tmix[v]);
            if (lower_limit >= upper_limit) return x0;
            return sample(x0);
        }


        double sample(double x0)
        {
            double d0;
            double lp0 = f(x0, d0);
            assert_finite(lp0);

            double slice_height = fastlog(
                    std::max<double>(random_uniform_01(rng), constants::zero_eps)) + lp0;

            double x_min_lp, x_max_lp;
            double x_min = find_slice_edge(x0, slice_height, lp0, d0, -1, &x_min_lp);
            double x_max = find_slice_edge(x0, slice_height, lp0, d0,  1, &x_max_lp);

            const double x_eps = 1e-6;
            double d;
            double x = (x_min + x_max) / 2;
            while (fabs(x_max - x_min) > x_eps) {
                x = x_min + (x_max - x_min) * random_uniform_01(rng);
                double lp = f(x, d);

                if (lp >= slice_height) break;
                else if (x > x0) x_max = x;
                else             x_min = x;
            }

            return x;
        }

        void prepare_sample_optimize(unsigned int u, unsigned int v)
        {
            this->u = u;
            this->v = v;
            double x0 = S.tmix[u] / (S.tmix[u] + S.tmix[v]);

            c = S.transcript_component[u];
            assert(c == S.transcript_component[v]);

            tgroup = S.transcript_tgroup[u];
            assert(tgroup == S.transcript_tgroup[v]);

            component_size = S.component_frag[c + 1] - S.component_frag[c];
            lp0 = dotlog(S.frag_counts[c], S.frag_probs[c], component_size) / M_LOG2E;

            wtgroupmix0 = 0.0;
            prior_lp0 = 0.0;
            if (S.use_priors) {
                BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
                    wtgroupmix0 += S.tmix[tid] / S.transcript_weights[tid];
                }

                double wtmixu = (S.tmix[u] / S.transcript_weights[u]) / wtgroupmix0;
                double wtmixv = (S.tmix[v] / S.transcript_weights[v]) / wtgroupmix0;

                prior_lp0 += splice_prior.f(S.hp.splice_mu[u], S.hp.splice_sigma[u],
                                            &wtmixu, 1);

                prior_lp0 += splice_prior.f(S.hp.splice_mu[v], S.hp.splice_sigma[v],
                                            &wtmixv, 1);

                prior_lp0 += log_weight_transform_gradient(x0);

                lp0 += prior_lp0;
                assert_finite(prior_lp0);
            }

            f0 = S.component_frag[c+1];
            f1 = S.component_frag[c];
            if (S.weight_matrix->rowlens[u] > 0) {
                f0 = S.weight_matrix->idxs[u][0];
                f1 = 1 + S.weight_matrix->idxs[u][S.weight_matrix->rowlens[u] - 1];
            }

            if (S.weight_matrix->rowlens[v] > 0) {
                f0 = std::min(f0, S.weight_matrix->idxs[v][0]);
                f1 = std::max(f1,  1 + S.weight_matrix->idxs[v][S.weight_matrix->rowlens[v] - 1]);
            }

            if (f0 > f1) f0 = f1 = S.component_frag[c];
            f0 -= S.component_frag[c];
            f1 -= S.component_frag[c];

            // Round f0 down to the nearst 32-byte boundry so we can use avx/sse.
            f0 &= 0xfffffff8;

            lpf01 = dotlog(S.frag_counts[c] + f0, S.frag_probs[c] + f0, f1 - f0) / M_LOG2E;
            assert_finite(lpf01);
        }


        double optimize(unsigned int u, unsigned int v)
        {
            prepare_sample_optimize(u, v);
            double x0 = S.tmix[u] / (S.tmix[u] + S.tmix[v]);
            if (lower_limit >= upper_limit) return x0;
            return optimize(x0);
        }


        double optimize(double x0)
        {
            double maxf;
            x0 = std::max<double>(std::min<double>(x0, upper_limit), lower_limit);

            nlopt_result result = nlopt_optimize(opt, &x0, &maxf);
            if (result < 0 && (result != NLOPT_FAILURE ||
                        !boost::math::isfinite(x0) || !boost::math::isfinite(maxf))) {
                Logger::warn("Optimization failed with code %d", (int) result);
            }
            assert_finite(x0);

            return x0;
        }

        double find_slice_edge(double x0, double slice_height,
                               double lp0, double d0, int direction,
                               double* last_lp)
        {
            const double lp_eps = 1e-3;
            const double d_eps  = 1e-1;
            const double x_eps  = 1e-6;

            double lp = lp0 - slice_height;
            double d = d0;
            double x = x0;
            double x_bound_lower, x_bound_upper;
            if (direction < 0) {
                x_bound_lower = lower_limit;
                x_bound_upper = x0;
            }
            else {
                x_bound_lower = x0;
                x_bound_upper = upper_limit;
            }

            while (fabs(lp) > lp_eps && fabs(x_bound_upper - x_bound_lower) > x_eps) {
                double x1 = x - lp / d;
                if (fabs(d) < d_eps || !boost::math::isfinite(x1)) {
                    x1 = (x_bound_lower + x_bound_upper) / 2;
                }

                // if we are very close to the boundry, and this iteration moves us past
                // the boundry, just give up.
                if (direction < 0 && fabs(x - lower_limit) <= x_eps && (x1 < x || lp > 0.0)) break;
                if (direction > 0 && fabs(x - upper_limit) <= x_eps && (x1 > x || lp > 0.0)) break;

                // if we are moving in the wrong direction (i.e. toward the other root),
                // use bisection to correct course.
                if (direction < 0) {
                    if (lp > 0) x_bound_upper = x;
                    else        x_bound_lower = x;
                }
                else {
                    if (lp > 0) x_bound_lower = x;
                    else        x_bound_upper = x;
                }

                bool bisect = x1 < x_bound_lower + x_eps || x1 > x_bound_upper - x_eps;

                // try using the gradient
                if (!bisect) {
                    x = x1;
                    lp = f(x, d) - slice_height;
                    bisect = !boost::math::isfinite(lp) || !boost::math::isfinite(d);
                }

                // resort to binary search if we seem not to be making progress
                if (bisect) {
                    size_t iteration_count = 0;
                    while (true) {
                        x = (x_bound_lower + x_bound_upper) / 2;
                        lp = f(x, d) - slice_height;

                        //if (boost::math::isinf(lp) || boost::math::isinf(d)) {
                        //if (!boost::math::isfinite(lp) || !boost::math::isfinite(d)) {
                        if (!boost::math::isfinite(lp)) {
                            if (direction < 0) x_bound_lower = x;
                            else               x_bound_upper = x;
                        }
                        else break;

                        if (++iteration_count > 50) {
                            Logger::abort("Slice sampler edge finding is not making progress.");
                        }
                    }
                }

                assert_finite(lp);
                //assert_finite(d);
            }

            assert_finite(x);
            *last_lp = lp;

            return x;
        }

    protected:
        // Since the prior is specified over weighted transcript abundance, we
        // have to account for the transform, or else the results can be a bit
        // biased.
        double log_weight_transform_gradient(double x)
        {

            double wu = S.transcript_weights[u];
            double wv = S.transcript_weights[v];
            double c = S.tmix[u] + S.tmix[v];
            double a = wtgroupmix0 - S.tmix[u] / wu - S.tmix[v] / wv;

            // derivative of: (x*c/u) / (x*c/u + (1-x)*c/v + a))
            double su =
                (c * wu * wv * (a * wv + c)) /
                sq(-a * wu * wv + c * wu * x - c * wu - c * wv * x);

            double log_su = fastlog(su);
            assert_finite(su);

            // derivative of: ((1-x)*c/v) / (x*c/u + (1-x)*c/v + a))
            double sv =
                - (c * wu * wv * (a * wu + c)) /
                sq(-a * wu * wv + c * wu * x - c * wu - c * wv * x);

            double log_sv = fastlog(-sv);
            assert_finite(sv);

            return log_su + log_sv;
        }


        // used for the jacobian adjustment
        double derivative_log_weight_transform_gradient(double x)
        {
            double wu = S.transcript_weights[u];
            double wv = S.transcript_weights[v];
            double c = S.tmix[u] + S.tmix[v];
            double a = wtgroupmix0 - S.tmix[u] / wu - S.tmix[v] / wv;

            double d = 2 * (2 * c * (wu - wv)) /
                (a * wu * wv + c * (wu * (-x) + wu + wv * x));

            assert_finite(d);

            return d;
        }


        double f(double x, double& d)
        {
            x = std::min<double>(upper_limit, std::max<double>(lower_limit, x));

            double tmixu = x * (S.tmix[u] + S.tmix[v]);
            double tmixv = (1 - x) * (S.tmix[u] + S.tmix[v]);

            acopy(S.frag_probs_prop[c] + f0, S.frag_probs[c] + f0, (f1 - f0) * sizeof(float));

            asxpy(S.frag_probs_prop[c],
                  S.weight_matrix->rows[u],
                  tmixu - S.tmix[u],
                  S.weight_matrix->idxs[u],
                  S.component_frag[c],
                  S.weight_matrix->rowlens[u]);

            asxpy(S.frag_probs_prop[c],
                  S.weight_matrix->rows[v],
                  tmixv - S.tmix[v],
                  S.weight_matrix->idxs[v],
                  S.component_frag[c],
                  S.weight_matrix->rowlens[v]);

            d = 0.0;
            d += asxtydsz(S.frag_counts[c], S.weight_matrix->rows[u], S.frag_probs_prop[c],
                          S.weight_matrix->idxs[u], S.component_frag[c],
                          S.weight_matrix->rowlens[u]);

            d -= asxtydsz(S.frag_counts[c], S.weight_matrix->rows[v], S.frag_probs_prop[c],
                          S.weight_matrix->idxs[v], S.component_frag[c],
                          S.weight_matrix->rowlens[v]);

            d *= tmixu + tmixv;

            double prior_lp_delta = 0.0;
            if (S.use_priors) {
                prior_lp_delta -= prior_lp0;

                double wtgroupmix = wtgroupmix0 +
                    (tmixu - S.tmix[u]) / S.transcript_weights[u] +
                    (tmixv - S.tmix[v]) / S.transcript_weights[v];

                double wtmixu = (tmixu / S.transcript_weights[u]) / wtgroupmix;
                double wtmixv = (tmixv / S.transcript_weights[v]) / wtgroupmix;

                prior_lp_delta += splice_prior.f(S.hp.splice_mu[u],
                                                 S.hp.splice_sigma[u],
                                                 &wtmixu, 1);

                prior_lp_delta += splice_prior.f(S.hp.splice_mu[v],
                                                 S.hp.splice_sigma[v],
                                                 &wtmixv, 1);

                // jacobian
                prior_lp_delta += log_weight_transform_gradient(x);

                d += splice_prior.df_dx(S.hp.splice_mu[u],
                                        S.hp.splice_sigma[u],
                                        &wtmixu, 1);

                d -= splice_prior.df_dx(S.hp.splice_mu[v],
                                        S.hp.splice_sigma[v],
                                        &wtmixv, 1);

                d += derivative_log_weight_transform_gradient(x);
            }

            return lp0 - lpf01 +
                   dotlog(S.frag_counts[c] + f0, S.frag_probs_prop[c] + f0, f1 - f0) / M_LOG2E +
                   prior_lp_delta;
        }


    private:
        nlopt_opt opt;

        NormalLogPdf splice_prior;

        unsigned int u, v, c, tgroup;
        double lp0, lpf01, prior_lp0;
        unsigned int component_size;
        double wtgroupmix0;
        unsigned int f0, f1;

        Sampler& S;

        // Bounds on the parameter being sampled over
        double lower_limit, upper_limit;

        rng_t rng;
        boost::random::uniform_01<double> random_uniform_01;

        friend double intertranscriptsampler_nlopt_objective(
                    unsigned int _n,const double* _x, double* _grad, void* data);
};



double intertranscriptsampler_nlopt_objective(unsigned int _n, const double* _x,
                                              double* _grad, void* data)
{
    UNUSED(_n);
    InterTranscriptSampler* sampler = reinterpret_cast<InterTranscriptSampler*>(data);
    if (_grad) return sampler->f(_x[0], _grad[0]);
    else {
        double d;
        return sampler->f(_x[0], d);
    }
}


double component_nlopt_objective(unsigned int, const double*, double*, void*);


/* The interface that MaxPostThread and MCMCThread implement. */
class AbundanceSamplerThread
{
    public:
        AbundanceSamplerThread(Sampler& S,
                        Queue<ComponentBlock>& q,
                        Queue<int>& notify_queue)
            : inter_tgroup_sampler(S)
            , inter_transcript_sampler(S)
            , S(S)
            , q(q)
            , notify_queue(notify_queue)
            , optimize_state(false)
            , thread(NULL)
        {
            component_opt = nlopt_create(NLOPT_LN_SBPLX, 1);
            nlopt_set_max_objective(component_opt,
                                    component_nlopt_objective,
                                    reinterpret_cast<void*>(this));
            double lower_bound = 1e-3;
            nlopt_set_lower_bounds(component_opt, &lower_bound);
            nlopt_set_ftol_abs(component_opt, 1e-2);
        }

        ~AbundanceSamplerThread()
        {
            if (thread) delete thread;

            nlopt_destroy(component_opt);
        }

        // Sample expression of component c
        void sample_component(unsigned int c);

        // Optimize expression of component c
        void optimize_component(unsigned int c);

        // Sample transcript mixtures within component c
        void sample_intra_component(unsigned int c);

        // Sample tgroup mixtures with a fixed component
        void sample_inter_tgroup(unsigned int c, unsigned int u, unsigned int v);

        // Sample transcript mixtures within a fixed tgroup
        void sample_intra_tgroup(unsigned int tgroup);

        // Sample relative abundance of two transcripts within the same tgroup
        void sample_inter_transcript(unsigned int u, unsigned int v);

        // If set to true, optimize rather than sample.
        void set_optimize(bool state)
        {
            optimize_state = state;
        }

        void run()
        {
            ComponentBlock block;
            while (true) {
                block = q.pop();
                if (block.is_end_of_queue()) break;

                this->rng = block.rng;

                for (unsigned int c = block.u; c < block.v; ++c) {
                    sample_intra_component(c);
                    if (optimize_state) {
                        optimize_component(c);
                    }
                    else {
                        sample_component(c);
                    }
                }

                notify_queue.push(1);
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&AbundanceSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }


    private:
        float find_component_slice_edge(unsigned int c, float x0, float slice_height, float step);

        float compute_component_probability(unsigned int c, float x);

        float recompute_intra_component_probability(unsigned int u, unsigned int v,
                                              float tmixu, float tmixv,
                                              unsigned int f0, unsigned int f1,
                                              float pf01, float p0,
                                              float bu, float bv,
                                              float tgroupmix_u,
                                              float tgroupmix_v,
                                              float *d);

        float transcript_slice_sample_search(float slice_height,
                                             unsigned int u, unsigned int v,
                                             float z0, float p0, bool left,
                                             double bu, double bv,
                                             double wtgroupmix_u,
                                             double wtgroupmix_v);

        NormalLogPdf cmix_prior;

        InterTgroupSampler inter_tgroup_sampler;
        InterTranscriptSampler inter_transcript_sampler;
        rng_t* rng;
        boost::random::uniform_01<double> random_uniform_01;
        boost::random::uniform_int_distribution<unsigned int> random_uniform_int;
        Sampler& S;
        Queue<ComponentBlock>& q;
        Queue<int>& notify_queue;
        bool optimize_state;
        boost::thread* thread;

        // currenty component being sampled or optimized over
        unsigned int c;

        nlopt_opt component_opt;
        friend double component_nlopt_objective(unsigned int, const double*, double*, void*);
};


double component_nlopt_objective(unsigned int _n, const double* _x,
                                 double* _grad, void* data)
{
    UNUSED(_n);
    UNUSED(_grad);
    AbundanceSamplerThread* sampler =
        reinterpret_cast<AbundanceSamplerThread*>(data);

    return sampler->compute_component_probability(sampler->c, _x[0]);
}


float AbundanceSamplerThread::find_component_slice_edge(unsigned int c,
                                                        float x0, float slice_height,
                                                        float step)
{
    float component_epsilon = S.component_num_transcripts[c] * constants::zero_eps;

    const float eps = 1e-3;
    float x, y;
    do {
        x = std::max<float>(component_epsilon, x0 + step);
        y = compute_component_probability(c, x);
        step *= 2;
    } while (y > slice_height && x > component_epsilon);
    step /= 2;

    // binary search to find the edge
    double a, b;
    if (step < 0.0) {
        a = std::max<float>(component_epsilon, x0 + step);
        b = x0;
    }
    else {
        a = x0;
        b = std::max<float>(component_epsilon, x0 + step);
    }

    while (fabs(b - a) / ((b+a)/2) > eps) {
        double mid = (a + b) / 2;
        double w = compute_component_probability(c, mid);

        if (step < 0) {
            if (w > slice_height) b = mid;
            else                  a = mid;
        }
        else {
            if (w > slice_height) a = mid;
            else                  b = mid;
        }
    }

    return (a + b) / 2;
}


float AbundanceSamplerThread::compute_component_probability(unsigned int c, float cmixc)
{
    float lp = gamma_lnpdf(S.frag_count_sums[c] +
                           constants::tmix_prior_prec,
                           1.0, cmixc);

    if (S.use_priors) {
        // prior
        BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
            double x = S.tgroupmix[tgroup] * cmixc;
            double logx = fastlog(x) - fastlog(S.tgroup_scaling[tgroup]);
            double mu = S.hp.tgroup_mu[tgroup];
            double sigma = S.hp.tgroup_sigma[tgroup];

            double prior_lp = cmix_prior.f(mu, sigma, &logx, 1);
            prior_lp -= fastlog(x); // jacobian
            lp += prior_lp;

            if (!boost::math::isfinite(lp)) {
                Logger::abort("non-finite log-likelihood encountered");
            }
        }
    }

    return lp;
}


void AbundanceSamplerThread::sample_intra_component(unsigned int c)
{
    if (S.component_num_transcripts[c] <= 1) return;

    double component_tmix = 0.0;
    BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
        S.tgroupmix[tgroup] = 0.0;

        BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
            S.tgroupmix[tgroup] += S.tmix[tid];
        }
        if (!boost::math::isfinite(S.tgroupmix[tgroup]) || S.tgroupmix[tgroup] <= 0.0) {
            Logger::abort("A non-positive tgroup mixture occurred.");
        }
        component_tmix += S.tgroupmix[tgroup];

        BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
            S.tgroup_tmix[tid] = S.tmix[tid] / S.tgroupmix[tgroup];
        }
    }

    acopy(S.frag_probs_prop[c], S.frag_probs[c],
          (S.component_frag[c + 1] - S.component_frag[c]) * sizeof(float));

    if (S.component_tgroups[c].size() > 1) {
        for (unsigned int i = 0; i < S.component_tgroups[c].size(); ++i) {
            random_uniform_int.param(boost::random::uniform_int_distribution<unsigned int>::param_type(
                        0, S.component_tgroups[c].size() - 1));
            unsigned int u = random_uniform_int(*rng);

            unsigned int v = u;
            while (v == u) {
                v = random_uniform_int(*rng);
            }

            unsigned int tgroupu = 0, tgroupv = 0;
            unsigned int j = 0;
            BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
                if (u == j) tgroupu = tgroup;
                if (v == j) tgroupv = tgroup;
                ++j;
            }

            sample_inter_tgroup(c, tgroupu, tgroupv);
        }
    }

    // sample transcript abundance within tgroups
    BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
        sample_intra_tgroup(tgroup);
    }
}



void AbundanceSamplerThread::sample_inter_tgroup(unsigned int c, unsigned int u, unsigned int v)
{
    double x = inter_tgroup_sampler.sample(*rng, c, u, v, optimize_state);

    double tgroupmix_u = x * (S.tgroupmix[u] + S.tgroupmix[v]);
    double tgroupmix_v = (1 - x) * (S.tgroupmix[u] + S.tgroupmix[v]);

    double u_delta = tgroupmix_u / S.tgroupmix[u];
    double v_delta = tgroupmix_v / S.tgroupmix[v];

    assert_finite(tgroupmix_u);
    assert_finite(tgroupmix_v);

    S.tgroupmix[u] = tgroupmix_u;
    S.tgroupmix[v] = tgroupmix_v;

    BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
        double tmix_tid = std::max<float>(constants::zero_eps, S.tmix[tid] * u_delta);
        asxpy(S.frag_probs[c],
              S.weight_matrix->rows[tid],
              tmix_tid - S.tmix[tid],
              S.weight_matrix->idxs[tid],
              S.component_frag[c],
              S.weight_matrix->rowlens[tid]);
        S.tmix[tid] = tmix_tid;
        assert_finite(S.tmix[tid]);
    }

    BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
        double tmix_tid = std::max<float>(constants::zero_eps, S.tmix[tid] * v_delta);
        asxpy(S.frag_probs[c],
              S.weight_matrix->rows[tid],
              tmix_tid - S.tmix[tid],
              S.weight_matrix->idxs[tid],
              S.component_frag[c],
              S.weight_matrix->rowlens[tid]);
        S.tmix[tid] = tmix_tid;
        assert_finite(S.tmix[tid]);
    }

    size_t component_size = S.component_frag[c+1] - S.component_frag[c];
    for (size_t i = 0; i < component_size; ++i) {
        // if fragment probabilities are just slightly negative, adjust to
        // avoid NaNs.
        if (S.frag_probs[c][i] < 0) S.frag_probs[c][i] = 0.0;

        if (S.frag_counts[c][i] > 0 &&
            (S.frag_probs[c][i] <= 0 || !boost::math::isfinite(S.frag_probs[c][i]))) {
            Logger::abort("Non-positive fragment probability while doing inter tgroup sampling.");
        }
    }
}


void AbundanceSamplerThread::sample_intra_tgroup(unsigned int tgroup)
{
    if (S.tgroup_tids[tgroup].size() > 1) {
        for (unsigned int i = 0; i < S.tgroup_tids[tgroup].size() - 1; ++i) {
            random_uniform_int.param(boost::random::uniform_int_distribution<unsigned int>::param_type(
                        0, S.tgroup_tids[tgroup].size() - 1));
            unsigned int u = random_uniform_int(*rng);

            unsigned int v = u;
            while (v == u) {
                v = random_uniform_int(*rng);
            }

            sample_inter_transcript(S.tgroup_tids[tgroup][u],
                                    S.tgroup_tids[tgroup][v]);
        }
    }
}


void AbundanceSamplerThread::sample_inter_transcript(unsigned int u, unsigned int v)
{
    double x;
    if (optimize_state)  {
        x = inter_transcript_sampler.optimize(u, v);
    }
    else {
        x = inter_transcript_sampler.sample(u, v);
    }

    float tmixu = x * (S.tmix[u] + S.tmix[v]);
    tmixu = std::max<double>(tmixu, constants::zero_eps);
    tmixu = std::min<double>(tmixu, 1.0 - constants::zero_eps);
    assert_finite(tmixu);

    float tmixv = (1 - x) * (S.tmix[u] + S.tmix[v]);
    tmixv = std::max<double>(tmixv, constants::zero_eps);
    tmixv = std::min<double>(tmixv, 1.0 - constants::zero_eps);
    assert_finite(tmixv);

    float tmix_delta_u = tmixu - S.tmix[u];
    float tmix_delta_v = tmixv - S.tmix[v];

    unsigned int c = S.transcript_component[u];
    assert(c == S.transcript_component[v]);

    unsigned int tgroup = S.transcript_tgroup[u];
    assert(tgroup = S.transcript_tgroup[v]);

    asxpy(S.frag_probs[c],
          S.weight_matrix->rows[u],
          tmix_delta_u,
          S.weight_matrix->idxs[u],
          S.component_frag[c],
          S.weight_matrix->rowlens[u]);

    asxpy(S.frag_probs[c],
          S.weight_matrix->rows[v],
          tmix_delta_v,
          S.weight_matrix->idxs[v],
          S.component_frag[c],
          S.weight_matrix->rowlens[v]);

    size_t component_size = S.component_frag[c+1] - S.component_frag[c];
    for (size_t i = 0; i < component_size; ++i) {
        // if fragment probabilities are just slightly negative, adjust to
        // avoid NaNs.
        if (S.frag_probs[c][i] < 0) S.frag_probs[c][i] = 0.0;

        if (S.frag_counts[c][i] > 0 &&
            (S.frag_probs[c][i] <= 0 || !boost::math::isfinite(S.frag_probs[c][i]))) {
            Logger::abort("Non-positive fragment probability while doing inter transcript sampling.");
        }
    }

    S.tgroup_tmix[u] = tmixu / S.tgroupmix[tgroup];
    S.tgroup_tmix[v] = tmixv / S.tgroupmix[tgroup];

    S.tmix[u] = tmixu;
    S.tmix[v] = tmixv;
}


/* I should write a little abount what I'm doing with cmix. Cmix lies on a high
 * dimensional simplex, but I'm sort of punting a bit and treating each cmix as
 * independent. To do so, I keep track of the factor by which it was scaled
 * *last* round to put it on a simplex and scale it by that. Super-weird I know,
 * but the simplex is thousands of dimensions so as an approximation it isn't so
 * bad.
 */
void AbundanceSamplerThread::sample_component(unsigned int c)
{
    BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
        S.tgroupmix[tgroup] = 0.0;
        BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
            S.tgroupmix[tgroup] += S.tmix[tid];
        }
    }

    this->c = c;

    float x0 = S.cmix[c];

    float lp0 = compute_component_probability(c, x0);
    float slice_height = lp0 +
        fastlog(std::max<double>(constants::zero_eps, random_uniform_01(*rng)));
    float step = 1.0;

    float x_min = find_component_slice_edge(c, x0, slice_height, -step);
    float x_max = find_component_slice_edge(c, x0, slice_height, +step);

    float x;
    while (true) {
         x = x_min + (x_max - x_min) * random_uniform_01(*rng);
         float lp = compute_component_probability(c, x);

         if (lp >= slice_height) break;
         else if (x > x0) x_max = x;
         else             x_min = x;

         if (fabs(x_max - x_min) < 1e-8) break;
    }

    S.cmix[c] = x;
}


void AbundanceSamplerThread::optimize_component(unsigned int c)
{
    BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
        S.tgroupmix[tgroup] = 0.0;
        BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
            S.tgroupmix[tgroup] += S.tmix[tid];
        }
    }

    this->c = c;
    double x0 = S.frag_count_sums[c] +
                constants::tmix_prior_prec * S.component_num_transcripts[c];

    double maxf;
    double component_epsilon = S.component_num_transcripts[c] * 1e-6;
    nlopt_set_lower_bounds(component_opt, &component_epsilon);
    nlopt_result result = nlopt_optimize(component_opt, &x0, &maxf);
    if (result < 0 && (result != NLOPT_FAILURE ||
                !boost::math::isfinite(x0) || !boost::math::isfinite(maxf))) {
        Logger::warn("Optimization failed with code %d", (int) result);
    }

    S.cmix[c] = x0;
}


class MultireadSamplerThread
{
    public:
        MultireadSamplerThread(Sampler& S,
                Queue<MultireadBlock>& q,
                Queue<int>& notify_queue)
            : hillclimb(false)
            , S(S)
            , q(q)
            , notify_queue(notify_queue)
            , thread(NULL)
    {
    }

        ~MultireadSamplerThread()
        {
            if (thread) delete thread;
        }

        void run()
        {
            MultireadBlock block;
            while (true) {
                block = q.pop();
                if (block.is_end_of_queue()) break;
                rng = block.rng;
                for (; block.u < block.v; ++block.u) {
                    float sumprob = 0.0;
                    unsigned int k = S.multiread_num_alignments[block.u];

                    for (unsigned int i = 0; i < k; ++i) {

                        unsigned int c = S.multiread_alignments[block.u][i].component;
                        unsigned int f = S.multiread_alignments[block.u][i].frag;
                        float align_pr = S.multiread_alignments[block.u][i].align_pr;
                        S.frag_counts[c][f - S.component_frag[c]] = 0;
                        sumprob += align_pr *
                                   S.cmix[c] *
                                   S.frag_probs[c][f - S.component_frag[c]];
                        //sumprob += align_pr;
                    }

                    //if (hillclimb) {
                    if (true) {
                        for (unsigned int i = 0; i < k; ++i) {
                            unsigned int c = S.multiread_alignments[block.u][i].component;
                            unsigned int f = S.multiread_alignments[block.u][i].frag;
                            float align_pr = S.multiread_alignments[block.u][i].align_pr;
                            S.frag_counts[c][f - S.component_frag[c]] =
                                (align_pr * S.cmix[c] * S.frag_probs[c][f - S.component_frag[c]]) / sumprob;
                        }
                    }
                    else {
                        float r = sumprob * random_uniform_01(*rng);
                        unsigned int i;
                        for (i = 0; i < k; ++i) {
                            unsigned int c = S.multiread_alignments[block.u][i].component;
                            unsigned int f = S.multiread_alignments[block.u][i].frag;
                            float align_pr = S.multiread_alignments[block.u][i].align_pr;

                            float p = align_pr *
                                      S.cmix[c] *
                                      S.frag_probs[c][f - S.component_frag[c]];

                            //float p = align_pr;

                            if (r <= p) {
                                break;
                            }
                            else {
                                r -= p;
                            }
                        }

                        i = std::min(i, k - 1);

                        unsigned int c = S.multiread_alignments[block.u][i].component;
                        unsigned int f = S.multiread_alignments[block.u][i].frag;

                        S.frag_counts[c][f - S.component_frag[c]] = 1;
                    }
                }

                notify_queue.push(1);
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&MultireadSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

        /* When true, optimize probability rather than sample. */
        bool hillclimb;

    private:
        Sampler& S;
        Queue<MultireadBlock>& q;
        Queue<int>& notify_queue;
        boost::thread* thread;
        rng_t* rng;
        boost::random::uniform_01<double> random_uniform_01;
};


Sampler::Sampler(unsigned int rng_seed,
                 const char* bam_fn, const char* fa_fn,
                 TranscriptSet& ts, FragmentModel& fm,
                 bool use_priors)
    : ts(ts)
    , fm(fm)
    , use_priors(use_priors)
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

    weight_matrix = new WeightMatrix(ts.size());
    transcript_weights = new double [ts.size()];
    transcript_gc = new double [ts.size()];
    TSVec<FragIdxCount> nz_frag_counts;
    TSVec<MultireadFrag> multiread_frags;
    tgroup_scaling.resize(ts.num_tgroups());

    std::vector<FragWeightEstimationThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragWeightEstimationThread(
                    fm,
                    read_indexer,
                    *weight_matrix,
                    nz_frag_counts,
                    multiread_frags,
                    transcript_weights,
                    transcript_gc,
                    q));
        threads.back()->start();
    }

    Logger::debug("Loci: %lu", (unsigned long) intervals.size());

    std::string task_name = std::string("Estimating fragment weights (") +
                            std::string(bam_fn) +
                            std::string(")");

    sam_scan(intervals, bam_fn, fa_fn, task_name.c_str());

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    /* Free a little space. */
    fm.multireads.clear();
    for (size_t i = 0; i < constants::num_threads; ++i) delete threads[i];

    size_t oldncol;
    unsigned int* idxmap = weight_matrix->compact(&oldncol);
    Logger::debug("Weight-matrix dimensions: %lu x %lu",
            (unsigned long) weight_matrix->nrow,
            (unsigned long) weight_matrix->ncol);

    /* Update frag_count indexes */
    for (TSVec<FragIdxCount>::iterator i = nz_frag_counts.begin();
            i != nz_frag_counts.end(); ++i) {
        i->first = i->first < oldncol ? idxmap[i->first] : UINT_MAX;
    }

    /* Update multiread fragment indexes */
    for (TSVec<MultireadFrag>::iterator i = multiread_frags.begin();
            i != multiread_frags.end(); ++i) {
        i->frag_idx = i->frag_idx < oldncol ? idxmap[i->frag_idx] : UINT_MAX;
    }

    delete [] idxmap;

    /* Find connected components. */
    unsigned int N = weight_matrix->ncol + weight_matrix->nrow;

    /* ds is disjoint set data structure, where ds[i] holds a pointer to item
     * i's parent, and ds[i] = i, when the node is the root. */
    unsigned int* ds = new unsigned int [N];
    for (size_t i = 0; i < N; ++i) ds[i] = i;

    for (WeightMatrixIterator entry(*weight_matrix);
            entry != WeightMatrixIterator(); ++entry) {
        disjset_union(ds, weight_matrix->ncol + entry->i, entry->j);
    }

    /* Trancription groups cannot be split across components. */
    tgroup_tids = ts.tgroup_tids();
    BOOST_FOREACH (std::vector<unsigned int>& tids, tgroup_tids) {
        if (tids.empty()) continue;
        BOOST_FOREACH (unsigned int tid, tids) {
            disjset_union(ds,
                    weight_matrix->ncol + tids[0],
                    weight_matrix->ncol + tid);
        }
    }

    for (size_t i = 0; i < N; ++i) disjset_find(ds, i);

    /* Label components sequentially. */
    {
        boost::unordered_map<unsigned int, unsigned int> component_label;
        for (size_t i = 0; i < N; ++i) {
            boost::unordered_map<unsigned int, unsigned int>::iterator c;
            c = component_label.find(ds[i]);
            if (c == component_label.end()) {
                component_label.insert(std::make_pair(ds[i], component_label.size()));
            }
        }

        for (size_t i = 0; i < N; ++i) ds[i] = component_label[ds[i]];
        num_components = component_label.size();
    }

    Logger::debug("Components: %lu", (unsigned long) num_components);

    /* Label transcript components */
    component_tgroups.resize(num_components);
    component_num_transcripts = new unsigned int [num_components];
    std::fill(component_num_transcripts,
            component_num_transcripts + num_components,
            0);
    transcript_component = new unsigned int [ts.size()];
    transcript_tgroup = new unsigned int [ts.size()];
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        unsigned int c = ds[weight_matrix->ncol + t->id];
        transcript_tgroup[t->id] = t->tgroup;
        transcript_component[t->id] = c;
        component_num_transcripts[c]++;
        component_tgroups[c].insert(t->tgroup);
    }

    component_transcripts = new unsigned int* [num_components];
    for (size_t i = 0; i < num_components; ++i) {
        component_transcripts[i] = new unsigned int [component_num_transcripts[i]];
    }

    std::fill(component_num_transcripts,
            component_num_transcripts + num_components,
            0);

    for (size_t i = 0; i < ts.size(); ++i) {
        unsigned int c = ds[weight_matrix->ncol + i];
        component_transcripts[c][component_num_transcripts[c]++] = i;
    }

    /* Now, reorder columns by component. */
    unsigned int* idxord = new unsigned int [weight_matrix->ncol];
    for (size_t i = 0; i < weight_matrix->ncol; ++i) idxord[i] = i;
    std::sort(idxord, idxord + weight_matrix->ncol, ComponentCmp(ds));
    idxmap = new unsigned int [weight_matrix->ncol];
    for (size_t i = 0; i < weight_matrix->ncol; ++i) idxmap[idxord[i]] = i;
    delete [] idxord;

    std::sort(ds, ds + weight_matrix->ncol);
    weight_matrix->reorder_columns(idxmap);

    /* Update frag_count indexes */
    for (TSVec<FragIdxCount>::iterator i = nz_frag_counts.begin();
            i != nz_frag_counts.end(); ++i) {
        i->first = i->first < oldncol ? idxmap[i->first] : UINT_MAX;
    }
    std::sort(nz_frag_counts.begin(), nz_frag_counts.end());

    /* Build fragment count array. */
    component_frag = new unsigned int [num_components + 1];
    frag_counts = new float* [num_components];
    std::fill(frag_counts, frag_counts + num_components, (float*) NULL);
    frag_probs = new float* [num_components];
    std::fill(frag_probs, frag_probs + num_components, (float*) NULL);
    frag_probs_prop = new float * [num_components];
    std::fill(frag_probs_prop, frag_probs_prop + num_components, (float*) NULL);
    TSVec<FragIdxCount>::iterator fc = nz_frag_counts.begin();
    for (size_t i = 0, j = 0; i < num_components; ++i) {
        component_frag[i] = j;
        assert(ds[j] <= i);
        if (ds[j] > i) continue;
        size_t k;
        for (k = j; k < weight_matrix->ncol && ds[j] == ds[k]; ++k) {}
        size_t component_size = k - j;

        if (component_size == 0) continue;

        frag_counts[i] =
            reinterpret_cast<float*>(aalloc(component_size * sizeof(float)));
        std::fill(frag_counts[i], frag_counts[i] + component_size, 1.0f);

        frag_probs[i] =
            reinterpret_cast<float*>(aalloc(component_size * sizeof(float)));
        frag_probs_prop[i] =
            reinterpret_cast<float*>(aalloc(component_size * sizeof(float)));

        while (fc != nz_frag_counts.end() && fc->first < k) {
            frag_counts[i][fc->first - j] = (float) fc->second;
            ++fc;
        }

        j = k;
    }
    component_frag[num_components] = weight_matrix->ncol;

    /* Update multiread fragment indexes. */
    for (TSVec<MultireadFrag>::iterator i = multiread_frags.begin();
            i != multiread_frags.end(); ++i) {
        i->frag_idx = idxmap[i->frag_idx];
    }
    std::sort(multiread_frags.begin(), multiread_frags.end());

    num_multireads = 0;
    unsigned int num_alignments = 0;
    for (unsigned int i = 0; i < multiread_frags.size(); ) {
        unsigned int j;
        for (j = i; j < multiread_frags.size(); ++j) {
            if (multiread_frags[i].multiread_num != multiread_frags[j].multiread_num) break;
        }

        if (j - i > 1) {
            ++num_multireads;
            num_alignments += j - i;
        }
        i = j;
    }

    multiread_num_alignments = new unsigned int [num_multireads];
    std::fill(multiread_num_alignments,
            multiread_num_alignments + num_multireads, 0);
    multiread_alignments = new MultireadAlignment* [num_multireads];
    multiread_alignment_pool = new MultireadAlignment [num_alignments];

    unsigned int multiread_num = 0, alignment_num = 0;
    for (unsigned int i = 0; i < multiread_frags.size(); ) {
        unsigned int j;
        for (j = i; j < multiread_frags.size(); ++j) {
            if (multiread_frags[i].multiread_num != multiread_frags[j].multiread_num) break;
        }

        if (j - i > 1) {
            multiread_alignments[multiread_num] =
                &multiread_alignment_pool[alignment_num];
            multiread_num_alignments[multiread_num] = j - i;

            for (; i < j; ++i) {
                unsigned int c = ds[multiread_frags[i].frag_idx];
                multiread_alignment_pool[alignment_num].component = c;
                multiread_alignment_pool[alignment_num].frag =
                    multiread_frags[i].frag_idx;
                multiread_alignment_pool[alignment_num].align_pr =
                    multiread_frags[i].align_pr;
                ++alignment_num;
            }
            ++multiread_num;
        }

        i = j;
    }

    delete [] ds;
    delete [] idxmap;

    frag_count_sums = new float [num_components];
    tmix            = new double [weight_matrix->nrow];
    tgroupmix       = new double [tgroup_tids.size()];
    tgroup_tmix     = new double [weight_matrix->nrow];
    cmix            = new double [num_components];

    expr.resize(weight_matrix->nrow);

    // make sure no tiny transcript weights exist
    for (size_t i = 0; i < weight_matrix->nrow; ++i) {
        transcript_weights[i] = std::max<double>(constants::min_transcript_weight, transcript_weights[i]);
    }

    // initialize hyperparameters
    hp.scale = 1.0;
    hp.tgroup_mu.resize(ts.num_tgroups(), 0.0);
    hp.tgroup_sigma.resize(ts.num_tgroups(), 0.0);
    hp.splice_mu.resize(ts.size());
    hp.splice_sigma.resize(ts.size());
    std::fill(hp.splice_mu.begin(), hp.splice_mu.end(), 0.5);
    std::fill(hp.splice_sigma.begin(), hp.splice_sigma.end(), 0.5);

    // allocate sample threads
    for (size_t i = 0; i < constants::num_threads; ++i) {
        multiread_threads.push_back(
            new MultireadSamplerThread(*this, multiread_queue, multiread_notify_queue));
        abundance_threads.push_back(
            new AbundanceSamplerThread(*this, component_queue, component_notify_queue));
    }

    // allocate and seed rngs
    size_t num_multiread_rngs =
        num_multireads / constants::sampler_multiread_block_size + 1;
    multiread_rng_pool.resize(num_multiread_rngs);
    BOOST_FOREACH (rng_t& rng, multiread_rng_pool) {
        rng.seed(rng_seed++);
    }

    size_t num_abundance_rngs =
        num_components / constants::sampler_component_block_size + 1;
    abundance_rng_pool.resize(num_abundance_rngs);
    BOOST_FOREACH (rng_t& rng, abundance_rng_pool) {
        rng.seed(rng_seed++);
    }
}


void Sampler::engage_priors()
{
    use_priors = true;
}


void Sampler::disengage_priors()
{
    use_priors = false;
}


Sampler::~Sampler()
{
    delete [] tmix;
    delete [] cmix;
    delete [] tgroupmix;
    delete [] tgroup_tmix;
    delete [] transcript_component;
    delete [] transcript_tgroup;
    for (size_t i = 0; i < num_components; ++i) {
        afree(frag_counts[i]);
        afree(frag_probs[i]);
        afree(frag_probs_prop[i]);
        delete [] component_transcripts[i];
    }
    delete [] component_transcripts;
    delete [] component_num_transcripts;
    delete [] multiread_alignments;
    delete [] multiread_alignment_pool;
    delete [] multiread_num_alignments;
    delete [] frag_counts;
    delete [] frag_count_sums;
    delete [] frag_probs;
    delete [] frag_probs_prop;
    delete [] component_frag;
    delete weight_matrix;
    delete [] transcript_weights;
    delete [] transcript_gc;

    for (std::vector<MultireadSamplerThread*>::iterator i = multiread_threads.begin();
         i != multiread_threads.end(); ++i) {
        delete *i;
    }

    for (std::vector<AbundanceSamplerThread*>::iterator i = abundance_threads.begin();
         i != abundance_threads.end(); ++i) {
        delete *i;
    }
}


void Sampler::start()
{
    std::fill(tgroup_scaling.begin(), tgroup_scaling.end(), 1.0);

    /* Initial mixtures */
    for (unsigned int i = 0; i < num_components; ++i) {
        unsigned int max_frags = 0;
        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            unsigned int u = component_transcripts[i][j];
            if (weight_matrix->rowlens[u] > max_frags) {
                max_frags = weight_matrix->rowlens[u];
            }
        }

        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            tmix[component_transcripts[i][j]] = 1.0 / component_num_transcripts[i];
        }
    }

    std::fill(cmix, cmix + num_components, 1.0 / (float) num_components);
    init_multireads();
    init_frag_probs();
    update_frag_count_sums();

    /* Initial cmix */
    for (unsigned int c = 0; c < num_components; ++c) {
        cmix[c] =
            constants::tmix_prior_prec +
            frag_count_sums[c];
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        multiread_threads[i]->start();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        abundance_threads[i]->start();
    }
}


void Sampler::stop()
{
    for (size_t i = 0; i < constants::num_threads; ++i) {
        multiread_queue.push(MultireadBlock());
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        multiread_threads[i]->join();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        component_queue.push(ComponentBlock());
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        abundance_threads[i]->join();
    }
}


void Sampler::sample()
{
    sample_abundance();

    // adjust for transcript weight
    double expr_total = 0.0;
    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        double tmixi = tmix[i];
        //double tmixi = tmix[i] < 1e-5 ? 1e-10 :
                         //(tmix[i] > 1.0 - 1e-5 ? 1.0 - 1e-10 : tmix[i]);
        expr[i] = cmix[transcript_component[i]] * (tmixi / transcript_weights[i]);
        expr_total += expr[i];
    }

    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        expr[i] /= expr_total;
    }

    // super-useful diagnostics
#if 0
    {
        FILE* out = fopen("gc-length-expr.tsv", "w");
        fprintf(out, "transcript_id\texons\tgc\tlength\tweight\tcount\texpr\n");
        for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
            fprintf(out, "%s\t%lu\t%0.4f\t%lu\t%e\t%lu\t%e\n",
                    t->transcript_id.get().c_str(),
                    (unsigned long) t->size(),
                    transcript_gc[t->id],
                    (unsigned long) t->exonic_length(),
                    transcript_weights[t->id],
                    (unsigned long) weight_matrix->rowlens[t->id],
                    expr[t->id]);
        }
        fclose(out);
    }
#endif

#if 0
    // super-useful diagnostics
    {
        FILE* out = fopen("gc-length-expr.tsv", "w");
        fprintf(out, "gc\tlength\tweight\tcount\texpr\n");
        for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
            fprintf(out, "%0.4f\t%lu\t%e\t%lu\t%e\n",
                    transcript_gc[t->id],
                    (unsigned long) t->exonic_length(),
                    transcript_weights[t->id],
                    (unsigned long) weight_matrix->rowlens[t->id],
                    expr[t->id]);
        }
        fclose(out);
    }
#endif
    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        expr[i] *= hp.scale;
    }

    std::fill(tgroup_scaling.begin(), tgroup_scaling.end(), 0.0);
    for (unsigned int tgroup = 0; tgroup < tgroup_tids.size(); ++tgroup) {
        double scaled_expr = 0.0;
        double unscaled_expr = 0.0;
        BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
            scaled_expr += expr[tid];
            unscaled_expr += cmix[transcript_component[tid]] * tmix[tid];
        }
        tgroup_scaling[tgroup] = unscaled_expr / scaled_expr;
    }
}


void Sampler::optimize()
{
    for (std::vector<AbundanceSamplerThread*>::iterator i = abundance_threads.begin();
         i != abundance_threads.end(); ++i) {
        (*i)->set_optimize(true);
    }

    sample();

    for (std::vector<AbundanceSamplerThread*>::iterator i = abundance_threads.begin();
         i != abundance_threads.end(); ++i) {
        (*i)->set_optimize(false);
    }
}


const std::vector<double>& Sampler::state() const
{
    return expr;
}


void Sampler::clear_multireads()
{
    for (size_t i = 0; i < num_multireads; ++i) {
        size_t k = multiread_num_alignments[i];
        for (size_t j = 0; j < k; ++j) {
            size_t c = multiread_alignments[i][j].component;
            size_t f = multiread_alignments[i][j].frag;
            frag_counts[c][f - component_frag[c]] = 0;
        }
    }
}


void Sampler::sample_multireads()
{
    unsigned int r, i;
    for (r = 0, i = 0; r + constants::sampler_multiread_block_size < num_multireads;
            r += constants::sampler_multiread_block_size, ++i) {
        multiread_queue.push(MultireadBlock(
                    r, r + constants::sampler_multiread_block_size,
                    multiread_rng_pool[i]));
    }
    if (r < num_multireads) {
        multiread_queue.push(MultireadBlock(r, num_multireads, multiread_rng_pool[i]));
    }

    for (r = 0, i = 0; r + constants::sampler_multiread_block_size < num_multireads;
            r += constants::sampler_multiread_block_size, ++i) {
        multiread_notify_queue.pop();
    }
    if (r < num_multireads) {
        multiread_notify_queue.pop();
    }

    update_frag_count_sums();
}


void Sampler::update_frag_count_sums()
{
    // Recompute total fragments in each component
    std::fill(frag_count_sums, frag_count_sums + num_components, 0.0f);
    total_frag_count = 0;
    for (size_t i = 0; i < num_components; ++i) {
        unsigned int component_size = component_frag[i + 1] - component_frag[i];
        for (unsigned int j = 0; j < component_size; ++j) {
            frag_count_sums[i] += frag_counts[i][j];
        }
        total_frag_count += frag_count_sums[i];
    }
}


void Sampler::sample_abundance()
{
    unsigned int c, i;
    for (c = 0, i = 0; c + constants::sampler_component_block_size < num_components;
            c += constants::sampler_component_block_size, ++i)
    {
        component_queue.push(ComponentBlock(
            c, c + constants::sampler_component_block_size,
            abundance_rng_pool[i]));
    }
    if (c < num_components) {
        component_queue.push(ComponentBlock(c, num_components, abundance_rng_pool[i]));
    }

    for (c = 0, i = 0; c + constants::sampler_component_block_size < num_components;
            c += constants::sampler_component_block_size, ++i)
    {
        component_notify_queue.pop();
    }
    if (c < num_components) {
        component_notify_queue.pop();
    }

    // check for numerical errors
    for (unsigned int i = 0; i < num_components; ++i) {
        for (unsigned int j = 0; j < component_frag[i + 1] - component_frag[i]; ++j ) {
            if (!boost::math::isfinite(frag_probs[i][j])) {
                Logger::warn("numerical error: %f", frag_probs[i][j]);
            }
        }
    }
}


void Sampler::init_multireads()
{
    for (unsigned int i = 0; i < num_multireads; ++i) {
        for (unsigned int j = 0; j < multiread_num_alignments[i]; ++j) {
            unsigned int c = multiread_alignments[i][j].component;
            unsigned int f = multiread_alignments[i][j].frag;
            frag_counts[c][f - component_frag[c]] = 0;
        }
    }
}


void Sampler::init_frag_probs()
{
    for (unsigned int i = 0; i < num_components; ++i) {
        unsigned component_size = component_frag[i + 1] - component_frag[i];
        std::fill(frag_probs[i], frag_probs[i] + component_size, 0.0);
    }

    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        unsigned int c = transcript_component[i];
        asxpy(frag_probs[c],
              weight_matrix->rows[i],
              tmix[i],
              weight_matrix->idxs[i],
              component_frag[c],
              weight_matrix->rowlens[i]);

    }

    // just checking
    for (unsigned int i = 0; i < num_components; ++i) {
        unsigned component_size = component_frag[i + 1] - component_frag[i];
        for (unsigned int j = 0; j < component_size; ++j) {
            if (!boost::math::isfinite(frag_probs[i][j]) || frag_probs[i][j] <= 0.0) {
                Logger::abort("Non-positive initial fragment probability.");
            }
        }
    }
}


unsigned long Sampler::num_frags() const
{
    return lround(total_frag_count);
}


unsigned long Sampler::num_alignments() const
{
    return weight_matrix->ncol;
}


