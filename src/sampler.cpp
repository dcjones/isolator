
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <climits>
#include <ctime>

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


/* A (fragment indx, count) pair. */
typedef std::pair<unsigned int, unsigned int> FragIdxCount;


/* A row compressed sparse matrix.
 */
class WeightMatrix
{
    public:
        WeightMatrix(unsigned int nrow, unsigned int max_ncol)
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

            align_prs = new float* [nrow];
            std::fill(align_prs, align_prs + nrow, (float*) NULL);

            align_pr_sum = new float [max_ncol];
            std::fill(align_pr_sum, align_pr_sum + max_ncol, 0.0f);

            this->max_ncol = max_ncol;

            compacted = false;
        }

        size_t memory_used() const
        {
            size_t nbytes = sizeof(float) * max_ncol +
                nrow * 2 * sizeof(unsigned int);

            for (size_t i = 0; i < nrow; ++i) {
                nbytes += reserved[i] * (2 * sizeof(float) + sizeof(unsigned int));
            }
            return nbytes;
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
                    free(align_prs[i]);
                }
            }
            delete [] rows;
            delete [] idxs;
            delete [] align_prs;
            delete [] align_pr_sum;
        }

        void push(unsigned int i, unsigned int j, float w, float align_pr)
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
                    newsize = reserved[i] + 1000;
                }

                rows[i] = realloc_or_die(rows[i], newsize);
                idxs[i] = realloc_or_die(idxs[i], newsize);
                align_prs[i] = realloc_or_die(align_prs[i], newsize);
                reserved[i] = newsize;
            }

            idxs[i][rowlens[i]] = j;
            rows[i][rowlens[i]] = w;
            align_prs[i][rowlens[i]] = align_pr;
            ++rowlens[i];
        };

        void increase_align_pr_sum(unsigned int j, float align_pr)
        {
            boost::lock_guard<boost::mutex> lock(
                    align_pr_sum_mutex[j % (sizeof(align_pr_sum_mutex) / sizeof(align_pr_sum_mutex[0]))]);

            align_pr_sum[j] += align_pr;
        }


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
            /* Compute ncol */
            ncol = 0;
            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    ncol = std::max<unsigned int>(ncol, idxs[i][j]);
                }
            }

            /* Normalize alignment probabilities. */
            for (unsigned int i = 0; i < nrow; ++i) {
                for (unsigned int j = 0; j < rowlens[i]; ++j) {
                    align_prs[i][j] /= align_pr_sum[idxs[i][j]];
                }
            }
            delete [] align_pr_sum;
            align_pr_sum = NULL;

            unsigned int total_supress_count = 0;

            /* Reallocate each row array */
            for (unsigned int i = 0; i < nrow; ++i) {
                unsigned int m = rowlens[i];
                if (m == 0) continue;

                unsigned int supressed_count = 0;
                for (unsigned int j = 0; j < m; ++j) {
                    if (align_prs[i][j] < constants::min_align_pr) {
                        ++supressed_count;
                        ++total_supress_count;
                    }
                }
                unsigned int newm = m - supressed_count;

                float* newrow = reinterpret_cast<float*>(
                        aalloc(newm * sizeof(float)));

                unsigned int* newidx = reinterpret_cast<unsigned int*>(
                        aalloc(newm * sizeof(unsigned int)));

                unsigned int k = 0;
                for (unsigned int j = 0; j < m; ++j) {
                    if (align_prs[i][j] >= constants::min_align_pr) {
                        newrow[k] = rows[i][j] * align_prs[i][j];
                        newidx[k] = idxs[i][j];
                        ++k;
                    }
                }

                free(rows[i]);
                rows[i] = newrow;

                free(idxs[i]);
                idxs[i] = newidx;

                reserved[i] = newm;
                rowlens[i] = newm;

                free(align_prs[i]);
                align_prs[i] = NULL;
            }

            Logger::debug("%lu low probability alignments supressed",
                         total_supress_count);

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

        void expand_row(float* out, unsigned int i, unsigned int off)
        {
            for (unsigned int j = 0; j < rowlens[i]; ++j) {
                out[idxs[i][j] - off] = rows[i][j];
            }
        }

        void expand_row(float* out, unsigned int i, unsigned int off, float w)
        {
            for (unsigned int j = 0; j < rowlens[i]; ++j) {
                out[idxs[i][j] - off] += w * rows[i][j];
            }
        }

        void expand_row(float* out, const unsigned int* idxs, unsigned int numidxs,
                        unsigned int i, float w = 1.0)
        {
            unsigned int j = 0, k = 0;
            while (k < numidxs) {
                if (j >= rowlens[i] || idxs[k] < this->idxs[i][j]) {
                    ++k;
                }
                else if (idxs[k] > this->idxs[i][j]) {
                    Logger::abort("Missing index in expand_row");
                }
                else {
                    out[k++] += w * rows[i][j++];
                }
            }
        }

        unsigned int idx_union(unsigned int* out, unsigned int u, unsigned int v)
        {
            return idx_union(out, idxs[u], rowlens[u], v);
        }

        unsigned int idx_union(unsigned int* out, const unsigned int* idxs_u,
                               unsigned int idxs_u_len, unsigned int v)
        {
            unsigned int ku = 0, kv = 0, k = 0;
            while (ku < idxs_u_len || kv < rowlens[v]) {
                if (ku >= idxs_u_len) {
                    out[k++] = idxs[v][kv++];
                }
                else if (kv >= rowlens[v]) {
                    out[k++] = idxs_u[ku++];
                }
                else if (idxs_u[ku] < idxs[v][kv]) {
                    out[k++] = idxs_u[ku++];
                }
                else if (idxs_u[ku] > idxs[v][kv]) {
                    out[k++] = idxs[v][kv++];
                }
                else {
                    out[k++] = idxs_u[ku++];
                    ++kv;
                }
            }
            return k;
        }

        unsigned int nrow, ncol, max_ncol;

        /* Number of entries in each row. */
        unsigned int* rowlens;

        /* Row data. */
        float** rows;

        /* Columns indexed for each row. */
        unsigned int** idxs;

        /* Alignment probabilities, used when building this matrix. */
        float** align_prs;

        /* align_pr normalization indexed by read id */
        float* align_pr_sum;

        /* mutexes for blocks of values in align_pr_sum */
        boost::mutex align_pr_sum_mutex[32];

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
            start = ts.min_start;
            end = ts.max_end;
        }

        void add_alignment(long idx, const bam1_t* b)
        {
            rs.add_alignment(idx, b);
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
                             std::set<std::string> excluded_seqs,
                             AlnIndex& alnindex,
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
              std::set<std::string> excluded_seqs,
              AlnIndex& alnindex,
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
    bool excluded_seq = false;
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
        if ((b->core.flag & BAM_FPAIRED) && (b->core.flag & BAM_FPROPER_PAIR) == 0) continue;

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort(
                    "The input SAM/BAM file must be sorted. "
                    "Please run: 'samtools sort'.");
        }

        if (b->core.tid != last_tid) {
            std::string seqname(bam_f->header->target_name[b->core.tid]);
            excluded_seq = excluded_seqs.find(seqname) != excluded_seqs.end();
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

        if (excluded_seq) continue;

        long idx = alnindex.get(bam1_qname(b));

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
                intervals[j]->add_alignment(idx, b);
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
                WeightMatrix& weight_matrix,
                double* transcript_weights,
                bool run_frag_correction,
                Queue<SamplerInitInterval*>& q)
            : weight_matrix(weight_matrix)
            , transcript_weights(transcript_weights)
            , fm(fm)
            , run_frag_correction(run_frag_correction)
            , q(q)
            , thread(NULL)
            , seqbias_size(0)
            , posbias_size(0)
        {
            seqbias[0] = seqbias[1] = NULL;
            posbias = NULL;
            tpbias = NULL;
            fragbias = NULL;
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
            afree(tpbias);
            afree(fragbias);
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
                                const Transcript& t,
                                bool run_frag_correction);

        /* Compute (a number proportional to) the probability of observing the
         * given fragment from the given transcript, assuming every transcript
         * is equally expressed. */
        float fragment_weight(
                const Transcript& t,
                const AlignmentPair& a);

        WeightMatrix& weight_matrix;
        double* transcript_weights;

        FragmentModel& fm;

        // True if fragmentation bias should be used
        bool run_frag_correction;

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
        float* tpbias;
        float* fragbias;
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
    std::vector<float> max_aln_prs(locus->rs.size());
    std::fill(max_aln_prs.begin(), max_aln_prs.end(), 0.0f);

    for (TranscriptSetLocus::iterator t = locus->ts.begin();
         t != locus->ts.end(); ++t) {

        pos_t L = constants::seqbias_left_pos, R = constants::seqbias_right_pos;
        if (locus->seq) {
            t->get_sequence(tseq0, *locus->seq, L, R);
            t->get_sequence(tseq1, *locus->seq, R, L);
            tseq1.revcomp();
        }

        transcript_sequence_bias(*locus, *t);
        transcript_weights[t->id] = transcript_weight(*locus, *t,
                                                      run_frag_correction);

        unsigned int i = 0;
        for (boost::unordered_map<long, AlignedRead*>::iterator r = locus->rs.rs.begin();
             r != locus->rs.rs.end(); ++r) {
            long idx = r->first;
            float max_w_align_pr = 0.0;
            float maxw = 0.0;
            float max_align_pr = 0.0;
            for (AlignedReadIterator a(*r->second); a != AlignedReadIterator(); ++a) {
                if ((a->mate1 && a->mate1->paired && !a->mate2) ||
                    (a->mate2 && a->mate2->paired && !a->mate1)) continue;

                float w = fragment_weight(*t, *a);
                float align_pr = 1.0;
                if (a->mate1) align_pr *= 1.0 - a->mate1->misaligned_pr;
                if (a->mate2) align_pr *= 1.0 - a->mate2->misaligned_pr;

                if (align_pr * w > max_w_align_pr) {
                    maxw = w;
                    max_align_pr = align_pr;
                    max_w_align_pr = align_pr * w;
                }
            }

            if (maxw > 0.0) {
                max_aln_prs[i] = std::max<float>(max_aln_prs[i], max_align_pr);
                weight_matrix.push(t->id, idx, maxw, max_align_pr);
            }
            ++i;
        }
    }

    unsigned int i = 0;
    for (boost::unordered_map<long, AlignedRead*>::iterator r = locus->rs.rs.begin();
         r != locus->rs.rs.end(); ++r, ++i) {
        long idx = r->first;
        weight_matrix.increase_align_pr_sum(idx, max_aln_prs[i]);
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
        return;
    }

    pos_t L = constants::seqbias_left_pos;

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
            posbias[pos - frag_len + 1] = fm.gcbias->get_bias(gc);
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
    const SamplerInitInterval& locus, const Transcript& t,
    bool run_frag_correction)
{
    pos_t trans_len = t.exonic_length();
    if ((size_t) trans_len + 1 > ws.size()) ws.resize(trans_len + 1);

    if ((size_t) trans_len > posbias_size) {
        afree(posbias);
        afree(tpbias);
        afree(fragbias);
        posbias = reinterpret_cast<float*>(aalloc(trans_len * sizeof(float)));
        tpbias = reinterpret_cast<float*>(aalloc(trans_len * sizeof(float)));
        fragbias = reinterpret_cast<float*>(aalloc(trans_len * sizeof(float)));
        posbias_size = trans_len;
    }

    // 3' bias
    {
        double omp = 1.0 - fm.tpbias->p;
        double c = 1.0;
        if (t.strand == strand_pos) {
            for (pos_t pos = trans_len - 1; pos >= 0; --pos) {
                c *= omp;
                tpbias[pos] = c;
            }
        } else {
            for (pos_t pos = 0; pos < trans_len; ++pos) {
                c *= omp;
                tpbias[pos] = c;
            }
        }
    }

    // fragmentation bias
    if (run_frag_correction) {
        fm.fragbias->get_bias(fragbias, trans_len);
        for (pos_t u = 0; u < trans_len; ++u) {
            fragbias[u] /= frag_len_c(trans_len - u);
        }
    }
    else {
        std::fill(fragbias, fragbias + trans_len, 1.0);
    }

    // combine fragbias and tpbias
    for (pos_t u = 0; u < trans_len; ++u) {
        tpbias[u] *= fragbias[u];
        if (!boost::math::isfinite(tpbias[u])) tpbias[u] = 0.0;
    }

    /* Set ws[k] to be the the number of fragmens of length k, weighted by
     * sequence bias. */
    memset(posbias, 1.0, trans_len * sizeof(float));
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = frag_len_p(frag_len);

        /* Don't bother considering sequence bias if the fragment length
         * probability is extremely small (i.e., so that it suffocates any
         * effect from bias.). */
        if (frag_len_pr < constants::transcript_len_min_frag_pr) {
            ws[frag_len] = (float) (trans_len - frag_len + 1);
            continue;
        }

        transcript_gc_bias(locus, t, frag_len);

        ws[frag_len] = 0.0;

        ws[frag_len] +=
            fm.strand_specificity *
            dot(&seqbias[0][0],
                &seqbias[1][frag_len - 1],
                posbias,
                tpbias,
                trans_len - frag_len + 1);

        ws[frag_len] +=
            (1.0 - fm.strand_specificity) *
            dot(&seqbias[1][0],
                &seqbias[0][frag_len - 1],
                posbias,
                tpbias,
                trans_len - frag_len + 1);
    }

    tw = 0.0;
    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = frag_len_p(frag_len);
        tw += frag_len_pr * ws[frag_len];
    }

    if (!boost::math::isfinite(tw) || tw <= constants::min_transcript_weight) {
        tw = constants::min_transcript_weight;
    }

    // fragment_weight expects this to be in genome orentation
    if (t.strand == strand_neg) std::reverse(tpbias, tpbias + trans_len);

    return tw;
}


float FragWeightEstimationThread::fragment_weight(const Transcript& t,
                                                  const AlignmentPair& a)
{
    pos_t frag_len = a.frag_len(t);
    pos_t L = constants::seqbias_left_pos;
    std::pair<pos_t, pos_t> clipping = a.flanks_soft_clipping();

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
    pos_t start_offset = 0, end_offset = 0;
    if (a.mate1 && a.mate2) {
        pos_t offset1, offset2;
        if (a.mate1->strand == strand_pos) {
            offset1 = t.get_offset(a.mate1->start, clipping.first, 0);
            offset2 = t.get_offset(a.mate2->end, 0, clipping.second);
            start_offset = offset1;
            end_offset = offset2;
        }
        else {
            offset1 = t.get_offset(a.mate1->end, 0, clipping.second);
            offset2 = t.get_offset(a.mate2->start, clipping.first, 0);
            start_offset = offset2;
            end_offset = offset1;
        }

        if (offset1 < 0 || offset1 >= tlen) return 0.0;
        w *= seqbias[a.mate1->strand == t.strand ? 0 : 1][offset1];

        if (offset2 < 0 || offset2 >= tlen) return 0.0;
        w *= seqbias[a.mate2->strand == t.strand ? 0 : 1][offset2];

        if (fm.gcbias) {
            float gc = tseq0.gc_count(
                    L + start_offset, L + end_offset) / (float) frag_len;
            w *= fm.gcbias->get_bias(gc / frag_len);
        }

        pos_t five_prime_end = a.mate1->strand == t.strand ? offset1 : offset2;
        w *= tpbias[five_prime_end];
    }
    else {
        if (a.mate1) {
            pos_t offset;
            if (a.mate1->strand == strand_pos) {
                offset = t.get_offset(a.mate1->start, clipping.first, 0);
                start_offset = offset;
                end_offset = std::min<pos_t>(tlen - 1, offset + frag_len - 1);
            }
            else {
                offset = t.get_offset(a.mate1->end, 0, clipping.second);
                start_offset = std::max<pos_t>(0, offset - frag_len + 1);
                end_offset = offset;
            }

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

            if (a.mate1->strand == t.strand) {
                w *= tpbias[offset];
            }
        }

        if (a.mate2) {
            pos_t offset;
            if (a.mate2->strand == strand_pos) {
                offset = t.get_offset(a.mate2->start, clipping.first, 0);
                start_offset = offset;
                end_offset = std::min<pos_t>(tlen - 1, offset + frag_len - 1);
            }
            else {
                offset = t.get_offset(a.mate2->end, 0, clipping.second);
                start_offset = std::max<pos_t>(0, offset - frag_len + 1);
                end_offset = offset;
            }

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

            if (a.mate2->strand == t.strand) {
                w *= tpbias[offset];
            }
        }
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

    if (frag_len_pr * w < constants::min_frag_weight) {
        return 0.0;
    }

    w = (w * frag_len_pr) / tw;

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


struct ComponentCmpRev
{
    ComponentCmpRev(unsigned int* ds)
        : ds(ds)
    {
    }

    bool operator () (unsigned int i, unsigned int j) {
        return ds[i] > ds[j];
    }

    unsigned int* ds;
};


class InterTgroupSampler : public Shredder
{
    public:
        InterTgroupSampler(Sampler& S)
            : Shredder(constants::zero_eps, 1.0 - constants::zero_eps, constants::zero_eps)
            , S(S)
        {
            size_t max_comp_size = 0;
            for (unsigned int c = 0; c < S.num_components; ++c) {
                max_comp_size = std::max<size_t>(max_comp_size,
                        S.component_frag[c + 1] - S.component_frag[c]);
            }

            row_u = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            row_v = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            frag_probs = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            frag_probs_prop = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            idxs = new unsigned int [max_comp_size];
            idxs_work = new unsigned int [max_comp_size];
        }

        ~InterTgroupSampler()
        {
            delete [] idxs;
            delete [] idxs_work;
            afree(frag_probs);
            afree(frag_probs_prop);
            afree(row_u);
            afree(row_v);
        }

        double sample(rng_t& rng, unsigned int c, unsigned int u, unsigned int v,
                      bool optimize_state)
        {
            double x0 = S.tgroupmix[u] / (S.tgroupmix[u] + S.tgroupmix[v]);

            if (upper_limit <= lower_limit) return x0;

            idxlen = 0;
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                std::swap(idxs, idxs_work);
                idxlen = S.weight_matrix->idx_union(idxs, idxs_work, idxlen, tid);
            }

            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                std::swap(idxs, idxs_work);
                idxlen = S.weight_matrix->idx_union(idxs, idxs_work, idxlen, tid);
            }

            uf0 = idxlen;
            uf1 = 0;
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                if (S.weight_matrix->rowlens[tid] > 0) {
                    unsigned int firstidx = S.weight_matrix->idxs[tid][0];
                    unsigned int lastidx = S.weight_matrix->idxs[tid][S.weight_matrix->rowlens[tid] - 1];

                    unsigned int tid_f0, tid_f1;
                    for (tid_f0 = 0; tid_f0 < idxlen && idxs[tid_f0] != firstidx; ++tid_f0);
                    for (tid_f1 = tid_f0; tid_f1 < idxlen && idxs[tid_f1] != lastidx; ++tid_f1);
                    ++tid_f1;

                    uf0 = std::min<unsigned int>(uf0, tid_f0);
                    uf1 = std::max<unsigned int>(uf1, tid_f1);
                }
            }
            uf0 &= 0xfffffff8;
            memset(row_u, 0, idxlen * sizeof(float));

            vf0 = idxlen;
            vf1 = 0;
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                if (S.weight_matrix->rowlens[tid] > 0) {
                    unsigned int firstidx = S.weight_matrix->idxs[tid][0];
                    unsigned int lastidx = S.weight_matrix->idxs[tid][S.weight_matrix->rowlens[tid] - 1];

                    unsigned int tid_f0, tid_f1;
                    for (tid_f0 = 0; tid_f0 < idxlen && idxs[tid_f0] != firstidx; ++tid_f0);
                    for (tid_f1 = tid_f0; tid_f1 < idxlen && idxs[tid_f1] != lastidx; ++tid_f1);
                    ++tid_f1;

                    vf0 = std::min<unsigned int>(vf0, tid_f0);
                    vf1 = std::max<unsigned int>(vf1, tid_f1);
                }
            }
            vf0 &= 0xfffffff8;

            this->x0 = x0;
            this->c = c;
            this->u = u;
            this->v = v;
            this->tgroupmix_uv = S.tgroupmix[u] + S.tgroupmix[v];

            memset(row_u, 0, idxlen * sizeof(float));
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[u]) {
                S.weight_matrix->expand_row(row_u, idxs + uf0, uf1 - uf0, tid,
                                            S.tgroup_tmix[tid]);
            }

            memset(row_v, 0, idxlen * sizeof(float));
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[v]) {
                S.weight_matrix->expand_row(row_v, idxs + vf0, vf1 - vf0, tid,
                                            S.tgroup_tmix[tid]);
            }

            for (unsigned int k = 0; k < idxlen; ++k) {
                frag_probs[k] = S.frag_probs[c][idxs[k] - S.component_frag[c]];
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

            acopy(frag_probs_prop, frag_probs, idxlen * sizeof(float));

            for (unsigned int k = 0; k < idxlen; ++k) {
                frag_probs_prop[k] = S.frag_probs[c][idxs[k] - S.component_frag[c]];
            }

            axpy(frag_probs_prop + uf0,
                 row_u, tgroupmix_u - S.tgroupmix[u], uf1 - uf0);

            axpy(frag_probs_prop + vf0,
                 row_v, tgroupmix_v - S.tgroupmix[v], vf1 - vf0);

            // gradient
            d = 0.0;
            d += sumdiv(row_u, frag_probs_prop + uf0, uf1 - uf0);
            d -= sumdiv(row_v, frag_probs_prop + vf0, vf1 - vf0);
            d *= tgroupmix_u + tgroupmix_v;

            // prior probability
            double prior_lp = 0.0;
            if (S.use_priors) {
                double xu = S.cmix[c] * tgroupmix_u;
                double logxu = fastlog(xu) - fastlog(S.tgroup_scaling[u]);
                double mu_u = S.hp.tgroup_mu[u];
                prior_lp += tgroup_prior.f(mu_u, S.hp.tgroup_sigma[u], logxu);
                prior_lp -= fastlog(x); // jacobian

                double xv = S.cmix[c] * tgroupmix_v;
                double logxv = fastlog(xv) - fastlog(S.tgroup_scaling[v]);
                double mu_v = S.hp.tgroup_mu[v];
                prior_lp += tgroup_prior.f(mu_v, S.hp.tgroup_sigma[v], logxv);
                prior_lp -= fastlog(1 - x); // jacobian

                // derivative of normal log-pdf
                d += (mu_u - logxu) / (sq(S.hp.tgroup_sigma[u]) * x);
                // XXX: Should this be +=, that makes more intuitive sense
                d -= (mu_v - logxv) / (sq(S.hp.tgroup_sigma[v]) * (1 - x));

                // derivative of jacobian
                d -= 1/x;
                d += 1/(1-x);
            }

            // optimization can fail if d is huge
            d = std::max<double>(std::min<double>(d, 1e4), -1e4);

            double lp = sumlog(frag_probs_prop, idxlen);
            return lp / M_LOG2E + prior_lp;
        }

    private:
        NormalLogPdf tgroup_prior;
        float *row_u, *row_v;
        unsigned int *idxs, *idxs_work;
        unsigned int idxlen;
        float *frag_probs, *frag_probs_prop;

        double x0;

        // log-probability of the fragments in [f0, f1]
        unsigned int c;
        unsigned int u;
        unsigned int v;
        unsigned int uf0, uf1,
                     vf0, vf1;
        double tgroupmix_uv;

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
            opt = nlopt_create(NLOPT_LD_SLSQP, 1);
            //opt = nlopt_create(NLOPT_LN_SBPLX, 1);
            nlopt_set_lower_bounds(opt, &lower_limit);
            nlopt_set_upper_bounds(opt, &upper_limit);
            nlopt_set_max_objective(opt, intertranscriptsampler_nlopt_objective,
                                    reinterpret_cast<void*>(this));
            nlopt_set_ftol_abs(opt, 1e-7);
            nlopt_set_maxeval(opt, 20);

            size_t max_comp_size = 0;
            for (unsigned int c = 0; c < S.num_components; ++c) {
                max_comp_size = std::max<size_t>(max_comp_size,
                        S.component_frag[c + 1] - S.component_frag[c]);
            }

            row_u = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            row_v = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            frag_probs = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            frag_probs_prop = reinterpret_cast<float*>(aalloc(max_comp_size * sizeof(float)));
            idxs = new unsigned int [max_comp_size];
        }

        ~InterTranscriptSampler()
        {
            nlopt_destroy(opt);
            delete [] idxs;
            afree(frag_probs);
            afree(frag_probs_prop);
            afree(row_u);
            afree(row_v);
        }

        double sample(unsigned int u, unsigned int v)
        {
            prepare_sample_optimize(u, v);
            double x0 = S.tmix[u] / (S.tmix[u] + S.tmix[v]);
            if (lower_limit >= upper_limit) return x0;
            double x = sample(x0);

            return x;
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

            c = S.transcript_component[u];
            assert(c == S.transcript_component[v]);

            tgroup = S.transcript_tgroup[u];
            assert(tgroup == S.transcript_tgroup[v]);

            wtgroupmix0 = 0.0;
            BOOST_FOREACH (unsigned int tid, S.tgroup_tids[tgroup]) {
                wtgroupmix0 += S.tmix[tid] / S.transcript_weights[tid];
            }

            idxlen = S.weight_matrix->idx_union(idxs, u, v);

            uf0 = uf1 = 0;
            if (S.weight_matrix->rowlens[u] > 0) {
                unsigned int firstidx = S.weight_matrix->idxs[u][0];
                unsigned int lastidx = S.weight_matrix->idxs[u][S.weight_matrix->rowlens[u] - 1];
                for (uf0 = 0; uf0 < idxlen && idxs[uf0] != firstidx; ++uf0);
                for (uf1 = uf0; uf1 < idxlen && idxs[uf1] != lastidx; ++uf1);
                ++uf1;
            }
            uf0 &= 0xfffffff8;
            memset(row_u, 0, idxlen * sizeof(float));
            S.weight_matrix->expand_row(row_u, idxs + uf0, uf1 - uf0, u);

            vf0 = vf1 = 0;
            if (S.weight_matrix->rowlens[v] > 0) {
                unsigned int firstidx = S.weight_matrix->idxs[v][0];
                unsigned int lastidx = S.weight_matrix->idxs[v][S.weight_matrix->rowlens[v] - 1];
                for (vf0 = 0; vf0 < idxlen && idxs[vf0] != firstidx; ++vf0);
                for (vf1 = vf0; vf1 < idxlen && idxs[vf1] != lastidx; ++vf1);
                ++vf1;
            }
            vf0 &= 0xfffffff8;
            memset(row_v, 0, idxlen * sizeof(float));
            S.weight_matrix->expand_row(row_v, idxs + vf0, vf1 - vf0, v);

            for (unsigned int k = 0; k < idxlen; ++k) {
                frag_probs[k] = S.frag_probs[c][idxs[k] - S.component_frag[c]];
            }
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
            double x = std::max<double>(std::min<double>(x0, upper_limit), lower_limit);

            nlopt_result result = nlopt_optimize(opt, &x, &maxf);
            if (result < 0 && result != NLOPT_ROUNDOFF_LIMITED) {
                Logger::warn("Optimization failed with code %d", (int) result);
            }
            if (!boost::math::isfinite(x)) {
                return x0;
            }
            else return x;
        }

        double find_slice_edge(double x0, double slice_height,
                               double lp0, double d0, int direction,
                               double* last_lp)
        {
            const double lp_eps = 1e-3;
            const double d_eps  = 1e-1;
            const double x_eps  = 1e-6;

            size_t newton_count = 0;
            double lp = lp0 - slice_height;
            double d = d0;
            double x = x0;
            double x_bound_lower, x_bound_upper;
            if (direction < 0) {
                x_bound_lower = lower_limit;
                x_bound_upper = x0;
                if (f(0.0, d) > slice_height) {
                    return 0.0;
                }
            }
            else {
                x_bound_lower = x0;
                x_bound_upper = upper_limit;
                if (f(1.0, d) > slice_height) {
                    return 1.0;
                }
            }

            int slice_edge_search_count = 0;

            while (fabs(lp) > lp_eps && fabs(x_bound_upper - x_bound_lower) > x_eps) {
                ++slice_edge_search_count;
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

                bool bisect = newton_count > constants::max_newton_steps ||
                    x1 < x_bound_lower + x_eps || x1 > x_bound_upper - x_eps;

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
                else ++newton_count;

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
        double weight_transform_gradient_u(double x)
        {
            double wu = S.transcript_weights[u];
            double wv = S.transcript_weights[v];
            double c = S.tmix[u] + S.tmix[v];
            double a = wtgroupmix0 - S.tmix[u] / wu - S.tmix[v] / wv;

            // derivative of: (x*c/u) / (x*c/u + (1-x)*c/v + a))
            return
                (c * wu * wv * (a * wv + c)) /
                sq(-a * wu * wv + c * wu * x - c * wu - c * wv * x);
        }

        double weight_transform_gradient_v(double x)
        {

            double wu = S.transcript_weights[u];
            double wv = S.transcript_weights[v];
            double c = S.tmix[u] + S.tmix[v];
            double a = wtgroupmix0 - S.tmix[u] / wu - S.tmix[v] / wv;

            // derivative of: ((1-x)*c/v) / (x*c/u + (1-x)*c/v + a))
            return
                - (c * wu * wv * (a * wu + c)) /
                sq(-a * wu * wv + c * wu * x - c * wu - c * wv * x);
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

            acopy(frag_probs_prop, frag_probs, idxlen * sizeof(float));

            axpy(frag_probs_prop + uf0, row_u, tmixu - S.tmix[u], uf1 - uf0);
            axpy(frag_probs_prop + vf0, row_v, tmixv - S.tmix[v], vf1 - vf0);

            d = 0.0;
            d += sumdiv(row_u, frag_probs_prop + uf0, uf1 - uf0);
            d -= sumdiv(row_v, frag_probs_prop + vf0, vf1 - vf0);

            d *= tmixu + tmixv;

            double prior_lp = 0.0;
            if (S.use_priors) {
                double wtgroupmix = wtgroupmix0 +
                    (tmixu - S.tmix[u]) / S.transcript_weights[u] +
                    (tmixv - S.tmix[v]) / S.transcript_weights[v];

                double wtmixu = (tmixu / S.transcript_weights[u]) / wtgroupmix;
                double wtmixv = (tmixv / S.transcript_weights[v]) / wtgroupmix;

                prior_lp += splice_prior.f(S.hp.splice_mu[u],
                                           S.hp.splice_sigma[u],
                                           wtmixu);

                prior_lp += splice_prior.f(S.hp.splice_mu[v],
                                           S.hp.splice_sigma[v],
                                           wtmixv);

                double su = weight_transform_gradient_u(x),
                       sv = weight_transform_gradient_v(x);

                // jacobian
                prior_lp += fastlog(su) + fastlog(-sv);

                d += splice_prior.df_dx(S.hp.splice_mu[u],
                                        S.hp.splice_sigma[u],
                                        wtmixu) * su;

                d += splice_prior.df_dx(S.hp.splice_mu[v],
                                        S.hp.splice_sigma[v],
                                        wtmixv) * sv;

                // derivative of the log jacobian
                d += derivative_log_weight_transform_gradient(x);
            }

            // avoid optimization failing due to extremely large
            // gradients
            d = std::max<double>(std::min<double>(d, 1e4), -1e4);

            double lp = sumlog(frag_probs_prop, idxlen);
            return lp / M_LOG2E + prior_lp;
        }


    private:
        nlopt_opt opt;

        NormalLogPdf splice_prior;

        unsigned int u, v, c, tgroup;
        double wtgroupmix0;

        unsigned int uf0, uf1,
                     vf0, vf1;

        // dense copies of the sparse weight_matrix rows
        unsigned int* idxs;
        unsigned int idxlen;
        float *frag_probs, *frag_probs_prop;
        float* row_u;
        float* row_v;

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

                for (unsigned int l = block.u; l < block.v; ++l) {
                    unsigned int c = S.ordered_components[l];
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

        double compute_component_probability(unsigned int c, float x);

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

        LogNormalLogPdf cmix_prior;

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

    const float eps = 1e-4;
    double x, y;
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
        b = x0 + step;
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


double AbundanceSamplerThread::compute_component_probability(unsigned int c, float cmixc)
{
    double lp = gamma_lnpdf((S.component_frag[c + 1] - S.component_frag[c]) +
                           constants::tmix_prior_prec,
                           1.0, cmixc);

    if (S.use_priors) {
        // prior
        BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
            double x = S.tgroupmix[tgroup] * cmixc;
            double scaledx = x / S.tgroup_scaling[tgroup];
            double mu = S.hp.tgroup_mu[tgroup];
            double sigma = S.hp.tgroup_sigma[tgroup];

            double prior_lp = cmix_prior.f(mu, sigma, &scaledx, 1);
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

    // sample transcript abundance within tgroups
    BOOST_FOREACH (unsigned int tgroup, S.component_tgroups[c]) {
        sample_intra_tgroup(tgroup);
    }

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

            // don't try to sample if there's very little room to move
            if ((S.tgroupmix[tgroupu] <= 2.0 * constants::zero_eps &&
                 S.tgroupmix[tgroupv] <= 2.0 * constants::zero_eps) ||
                (S.tgroupmix[tgroupu] >= 1.0 - 2.0 * constants::zero_eps &&
                 S.tgroupmix[tgroupv] >= 1.0 - 2.0 * constants::zero_eps)) continue;

            sample_inter_tgroup(c, tgroupu, tgroupv);
        }
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
        double tmix_tid =
            std::min<float>(
                    1.0 - constants::zero_eps,
                    std::max<float>(constants::zero_eps, S.tmix[tid] * u_delta));
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
        double tmix_tid =
            std::min<float>(
                    1.0 - constants::zero_eps,
                    std::max<float>(constants::zero_eps, S.tmix[tid] * v_delta));
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

            unsigned int tidu = S.tgroup_tids[tgroup][u],
                         tidv = S.tgroup_tids[tgroup][v];

            // don't try to sample if there's very little room to move
            if ((S.tmix[tidu] <= 2.0 * constants::zero_eps &&
                 S.tmix[tidv] <= 2.0 * constants::zero_eps) ||
                (S.tmix[tidu] >= 1.0 - 2.0 * constants::zero_eps &&
                 S.tmix[tidv] >= 1.0 - 2.0 * constants::zero_eps)) continue;

            sample_inter_transcript(tidu, tidv);
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

    double lp0 = compute_component_probability(c, x0);
    double slice_height = lp0 +
        fastlog(std::max<double>(constants::zero_eps, random_uniform_01(*rng)));
    double step = 1.0;

    double x_min = find_component_slice_edge(c, x0, slice_height, -step);
    double x_max = find_component_slice_edge(c, x0, slice_height, +step);

    double x;
    while (true) {
         x = x_min + (x_max - x_min) * random_uniform_01(*rng);
         double lp = compute_component_probability(c, x);

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
    double x0 = (S.component_frag[c + 1] - S.component_frag[c]) +
                constants::tmix_prior_prec;
    double x = x0;

    double maxf;
    double component_epsilon = S.component_num_transcripts[c] * 1e-6;
    nlopt_set_lower_bounds(component_opt, &component_epsilon);
    nlopt_result result = nlopt_optimize(component_opt, &x, &maxf);
    if (result < 0 && result != NLOPT_ROUNDOFF_LIMITED) {
        Logger::warn("Optimization failed with code %d", (int) result);
    }

    S.cmix[c] = x;
}


Sampler::Sampler(unsigned int rng_seed,
                 const char* bam_fn, const char* fa_fn,
                 std::set<std::string> excluded_seqs,
                 TranscriptSet& ts, FragmentModel& fm,
                 bool run_frag_correction, bool use_priors)
    : ts(ts)
    , fm(fm)
    , run_frag_correction(run_frag_correction)
    , use_priors(use_priors)
{
    /* Producer/consumer queue of intervals containing indexed reads to be
     * processed. */
    Queue<SamplerInitInterval*> q(100);

    /* Collect intervals */
    std::vector<SamplerInitInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    weight_matrix = new WeightMatrix(ts.size(), fm.alnindex.size());
    transcript_weights = new double [ts.size()];
    tgroup_scaling.resize(ts.num_tgroups());

    std::vector<FragWeightEstimationThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragWeightEstimationThread(
                    fm,
                    *weight_matrix,
                    transcript_weights,
                    run_frag_correction,
                    q));
        threads.back()->start();
    }

    Logger::debug("Loci: %lu", (unsigned long) intervals.size());

    std::string task_name = std::string("Estimating fragment weights (") +
                            std::string(bam_fn) +
                            std::string(")");

    sam_scan(intervals, bam_fn, fa_fn, excluded_seqs, fm.alnindex, task_name.c_str());

    Logger::debug("weight matrix is %0.2fMB before compact",
                  (double) weight_matrix->memory_used() / 1e6);

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    /* Free a little space. */
    fm.alnindex.clear();
    for (size_t i = 0; i < constants::num_threads; ++i) delete threads[i];

    size_t oldncol;
    unsigned int* idxmap = weight_matrix->compact(&oldncol);
    Logger::debug("Weight-matrix dimensions: %lu x %lu",
            (unsigned long) weight_matrix->nrow,
            (unsigned long) weight_matrix->ncol);
    delete [] idxmap;

    Logger::debug("weight matrix is %0.2fMB after compact",
                  (double) weight_matrix->memory_used() / 1e6);

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

    /* Now, reorder columns by transcript and component. */

    unsigned int* idxord = new unsigned int [weight_matrix->ncol];
    for (size_t i = 0; i < weight_matrix->ncol; ++i) idxord[i] = i;

    std::stable_sort(idxord, idxord + weight_matrix->ncol, ComponentCmp(ds));

    idxmap = new unsigned int [weight_matrix->ncol];
    for (size_t i = 0; i < weight_matrix->ncol; ++i) idxmap[idxord[i]] = i;
    delete [] idxord;

    std::sort(ds, ds + weight_matrix->ncol);
    weight_matrix->reorder_columns(idxmap);

    /* Build fragment count array. */
    component_frag = new unsigned int [num_components + 1];
    frag_probs = new float* [num_components];
    std::fill(frag_probs, frag_probs + num_components, (float*) NULL);
    for (size_t i = 0, j = 0; i < num_components; ++i) {
        component_frag[i] = j;
        assert(ds[j] <= i);
        if (ds[j] > i) continue;
        size_t k;
        for (k = j; k < weight_matrix->ncol && ds[j] == ds[k]; ++k) {}
        size_t component_size = k - j;

        if (component_size == 0) continue;

        frag_probs[i] =
            reinterpret_cast<float*>(aalloc(component_size * sizeof(float)));

        j = k;
    }
    component_frag[num_components] = weight_matrix->ncol;

    delete [] ds;
    delete [] idxmap;

    tmix            = new double [weight_matrix->nrow];
    tgroupmix       = new double [tgroup_tids.size()];
    tgroup_tmix     = new double [weight_matrix->nrow];
    cmix            = new double [num_components];

    expr.resize(weight_matrix->nrow);

    // make sure no tiny transcript weights exist
    for (size_t i = 0; i < weight_matrix->nrow; ++i) {
        transcript_weights[i] = std::max<double>(constants::min_transcript_weight, transcript_weights[i]);
    }

    // figure out the best order to process components
    ordered_components.resize(num_components);
    for (unsigned int i = 0; i < num_components; ++i) {
        ordered_components[i] = i;
    }

    unsigned int* component_num_frags = new unsigned int [num_components];
    for (unsigned int c = 0; c < num_components; ++c) {
        component_num_frags[c] = component_frag[c+1] - component_frag[c];
    }

    std::sort(ordered_components.begin(), ordered_components.end(),
              ComponentCmpRev(component_num_frags));

    delete [] component_num_frags;

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
        abundance_threads.push_back(
            new AbundanceSamplerThread(*this, component_queue, component_notify_queue));
    }

    // allocate and seed rngs
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
        afree(frag_probs[i]);
        delete [] component_transcripts[i];
    }
    delete [] component_transcripts;
    delete [] component_num_transcripts;
    delete [] frag_probs;
    delete [] component_frag;
    delete weight_matrix;
    delete [] transcript_weights;

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
    init_frag_probs();

    /* Initial cmix */
    for (unsigned int c = 0; c < num_components; ++c) {
        cmix[c] =
            constants::tmix_prior_prec +
            (component_frag[c + 1] - component_frag[c]);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        abundance_threads[i]->start();
    }
}


void Sampler::stop()
{
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

    std::vector<double> cmix_weight(num_components);
    for (unsigned int i = 0; i < num_components; ++i) {
        cmix_weight[i] = 0.0;
        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            cmix_weight[i] += tmix[component_transcripts[i][j]] * transcript_weights[component_transcripts[i][j]];
        }
    }

    // adjust for transcript weight
    double expr_total = 0.0;
    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        double tmixi = tmix[i];
        //expr[i] = tmixi * cmix[transcript_component[i]] / cmix_weight[transcript_component[i]];
        expr[i] = tmixi * cmix[transcript_component[i]] / transcript_weights[i];
        expr_total += expr[i];
    }

    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        expr[i] /= expr_total;
    }

    // super-useful diagnostics
#if 0
    {
        FILE* out = fopen("gc-length-expr.tsv", "w");
        fprintf(out, "transcript_id\tstrand\texons\tlength\tweight\tcount\texpr\tcomp_size\tnum_tgroups\n");
        for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
            fprintf(out, "%s\t%d\t%lu\t%lu\t%e\t%lu\t%e\t%lu\t%lu\n",
                    t->transcript_id.get().c_str(),
                    (int) t->strand,
                    (unsigned long) t->size(),
                    (unsigned long) t->exonic_length(),
                    transcript_weights[t->id],
                    (unsigned long) weight_matrix->rowlens[t->id],
                    expr[t->id],
                    (unsigned long) component_num_transcripts[transcript_component[t->id]],
                    (unsigned long) component_tgroups[transcript_component[t->id]].size());
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
    return weight_matrix->ncol;
}


unsigned long Sampler::num_alignments() const
{
    return weight_matrix->ncol;
}


