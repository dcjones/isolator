
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <nlopt.h>

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
        unsigned int* compact()
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
            while (i < owner.nrow && k >= owner.rowlens[i]) {
                ++i;
                k = 0;
            }

            if (i < owner.nrow) {
                entry.i = i;
                entry.j = owner.idxs[i][k];
                entry.w = owner.rows[i][k];
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
                TSVec<FragIdxCount>& frag_counts,
                float* transcript_weights,
                Queue<SamplerInitInterval*>& q)
            : weight_matrix(weight_matrix)
            , frag_counts(frag_counts)
            , transcript_weights(transcript_weights)
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
        TSVec<FragIdxCount>& frag_counts;
        float* transcript_weights;

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

    for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
        /* fragment with count == 1 are implicit */
        if (f->second.second > 1) {
            frag_counts.push_back(f->second);
        }
    }

    std::vector<MultireadEntry> multiread_entries;
    for (TranscriptSetLocus::iterator t = locus->ts.begin();
            t != locus->ts.end(); ++t) {
        transcript_sequence_bias(*locus, *t);
        transcript_weights[t->id] = std::max(constants::min_transcript_weight,
                                             transcript_weight(*t));

        for (MultireadSet::iterator r = multiread_set.begin();
             r != multiread_set.end(); ++r) {
            for (AlignedReadIterator a(*r->second); a != AlignedReadIterator(); ++a) {
                float w = fragment_weight(*t, *a);
                if (w > constants::min_frag_weight) {
                    /* TODO: estimate alignment probability */
                    multiread_entries.push_back(
                            MultireadEntry(r->first, t->id, w, 1.0));
                }
            }
        }

        for (Frags::iterator f = frag_idx.begin(); f != frag_idx.end(); ++f) {
            float w = fragment_weight(*t, f->first);
            if (w > constants::min_frag_weight) {
                weight_matrix.push(t->id, f->second.first, w);
            }
        }
    }

    /* Process multiread entries into weight matrix entries */
    std::sort(multiread_entries.begin(), multiread_entries.end());
    for (std::vector<MultireadEntry>::iterator k, i = multiread_entries.begin();
            i != multiread_entries.end(); i = k) {
        k = i + 1;
        std::vector<MultireadEntry>::iterator j;
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
                if (w > constants::min_frag_weight) {
                    weight_matrix.push(j0->transcript_idx, j0->multiread_num, w);
                }
                j0 = j;
                w = 0.0;
            }

            w += j->align_pr / sum_align_pr * j->frag_weight;
        }
        if (w > constants::min_frag_weight) {
            weight_matrix.push(j0->transcript_idx, j0->multiread_num, w);
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

    weight_matrix = new WeightMatrix(ts.size());
    transcript_weights = new float [ts.size()];
    TSVec<FragIdxCount> nz_frag_counts;

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(
                    fm,
                    read_indexer,
                    *weight_matrix,
                    nz_frag_counts,
                    transcript_weights,
                    q));
        threads.back()->start();
    }

    Logger::info("Loci: %lu", (unsigned long) intervals.size());

    sam_scan(intervals, bam_fn, fa_fn, "Estimating fragment weights");

    /* Free a little space. */
    fm.multireads.clear();

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    unsigned int* idxmap = weight_matrix->compact();
    Logger::info("Weight-matrix dimensions: %lu x %lu",
            (unsigned long) weight_matrix->nrow,
            (unsigned long) weight_matrix->ncol);

    /* Update frag_count indexes */
    for (TSVec<FragIdxCount>::iterator i = nz_frag_counts.begin();
            i != nz_frag_counts.end(); ++i) {
        i->first = idxmap[i->first];
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

    Logger::info("Components: %lu", (unsigned long) num_components);

    /* Label transcript components */
    component_num_transcripts = new unsigned int [num_components];
    std::fill(component_num_transcripts,
              component_num_transcripts + num_components,
              0);
    transcript_component = new unsigned int [ts.size()];
    for (size_t i = 0; i < ts.size(); ++i) {
        unsigned int c = ds[weight_matrix->ncol + i];
        transcript_component[i] = c;
        component_num_transcripts[c]++;
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
        i->first = idxmap[i->first];
    }
    std::sort(nz_frag_counts.begin(), nz_frag_counts.end());

    /* Build fragment count array. */
    component_frag = new unsigned int [num_components + 1];
    frag_counts = new float* [num_components];
    std::fill(frag_counts, frag_counts + num_components, (float*) NULL);
    frag_probs = new float* [num_components];
    std::fill(frag_probs, frag_probs + num_components, (float*) NULL);
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

        while (fc != nz_frag_counts.end() && fc->first < k) {
            frag_counts[i][fc->first - j] = (float) fc->second;
            ++fc;
        }

        j = k;
    }
    component_frag[num_components] = weight_matrix->ncol;

    frag_count_sums = new float [num_components];
    std::fill(frag_count_sums, frag_count_sums + num_components, 0.0f);
    for (size_t i = 0; i < num_components; ++i) {
        unsigned int component_size = component_frag[i + 1] - component_frag[i];
        for (unsigned int j = 0; j < component_size; ++j) {
            frag_count_sums[i] += frag_counts[i][j];
        }
    }

    delete [] ds;
    delete [] idxmap;

    tmix = new double [weight_matrix->nrow];
    cmix = new double [num_components];
}


Sampler::~Sampler()
{
    delete [] tmix;
    delete [] cmix;
    delete [] transcript_component;
    for (size_t i = 0; i < num_components; ++i) {
        afree(frag_counts[i]);
        afree(frag_probs[i]);
        delete [] component_transcripts[i];
    }
    delete [] component_transcripts;
    delete [] component_num_transcripts;
    delete [] frag_counts;
    delete [] frag_count_sums;
    delete [] frag_probs;
    delete [] component_frag;
    delete weight_matrix;
    delete [] transcript_weights;
}


struct ComponentSubset
{
    ComponentSubset()
        : cs(NULL)
        , len(0)
    {
    }

    ComponentSubset(unsigned int* cs, unsigned int len)
        : cs(cs)
        , len(len)
    {
    }

    /* A sequence of components */
    unsigned int* cs;

    /* Number of components in the sequence */
    unsigned int len;

    /* Is this a special value marking the end of the queue. */
    bool is_end_of_queue() const
    {
        return len == 0;
    }
};



/* The interface that MaxPostThread and MCMCThread implement. */
class InferenceThread
{
    public:
        InferenceThread(Sampler& S,
                       Queue<ComponentSubset>& q)
            : S(S)
            , q(q)
            , thread(NULL)
        {
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            unsigned long seed = reinterpret_cast<unsigned long>(this) *
                                 (unsigned long) time(NULL);
            gsl_rng_set(rng, seed);
        }

        ~InferenceThread()
        {
            gsl_rng_free(rng);
            if (thread) delete thread;
        }

        /* Inference of relative expression of two transcripts a and b.. */
        virtual void run_inter_transcript(unsigned int a, unsigned int b) = 0;

        /* Inference relative expression of two components a and b. */
        virtual void run_inter_component(unsigned int a, unsigned int b) = 0;

        void run_intra_component(unsigned int c)
        {
            if (S.component_num_transcripts[c] == 0) return;

            gsl_ran_shuffle(rng, S.component_transcripts[c], 
                            S.component_num_transcripts[c],
                            sizeof(unsigned int));
            for (unsigned int i = 0; i < S.component_num_transcripts[c] - 1; ++i) {
                run_inter_transcript(S.component_transcripts[c][i],
                                     S.component_transcripts[c][i + 1]);
            }
        }

        void run()
        {
            ComponentSubset s;
            while (true) {
                s = q.pop();
                if (s.is_end_of_queue()) break;

                for (unsigned int i = 0; i < s.len; ++i) {
                    run_intra_component(s.cs[i]);
                }

                for (unsigned int i = 0; i < s.len - 1; ++i) {
                    run_inter_component(s.cs[i], s.cs[i + 1]);
                }
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&InferenceThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }


    protected:
        gsl_rng* rng;
        Sampler& S;

    private:
        Queue<ComponentSubset>& q;
        boost::thread* thread;
};


/* The parameter passed to nlopt objective functions. */
struct PairwiseObjFuncParam
{
    Sampler* S;
    unsigned int u, v;
};


double transcript_pair_objf(unsigned int n, const double* theta,
                            double* grad, void* param)
{
    assert(n == 1);
    assert(finite(*theta));
    UNUSED(n);

    PairwiseObjFuncParam* objfuncparam =
        reinterpret_cast<PairwiseObjFuncParam*>(param);
    unsigned int u = objfuncparam->u;
    unsigned int v = objfuncparam->v;
    Sampler& S = *objfuncparam->S;

    unsigned int c = S.transcript_component[u];
    assert(c == S.transcript_component[v]);

    double z = S.tmix[u] + S.tmix[v];
    S.tmix[u] = *theta * z;
    S.tmix[v] = (1.0 - *theta) * z;

    size_t component_size = S.component_frag[c + 1] - S.component_frag[c];
    std::fill(S.frag_probs[c], S.frag_probs[c] + component_size, 0.0);

    /* TODO: avoid totally recomputing fragment probabilities by subtracting the
     * previous value and adding the new one. */
    for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
        unsigned int t = S.component_transcripts[c][i];
        asxpy(S.frag_probs[c],
              S.weight_matrix->rows[t],
              S.tmix[t],
              S.weight_matrix->idxs[t],
              S.component_frag[c],
              S.weight_matrix->rowlens[t]);
    }

    double p = dotlog(S.frag_counts[c], S.frag_probs[c], component_size);
    assert(finite(p));

    if (grad != NULL) {
        *grad = 0.0;

        for (unsigned int j = 0; j < S.weight_matrix->rowlens[u]; ++j) {
            float denom = 
                S.frag_probs[c][
                    S.weight_matrix->idxs[u][j] - S.component_frag[c]];
            denom *= M_LN2;
            *grad += S.weight_matrix->rows[u][j] / denom;
            assert(finite(*grad));
        }

        for (unsigned int j = 0; j < S.weight_matrix->rowlens[v]; ++j) {
            float denom = 
                S.frag_probs[c][
                    S.weight_matrix->idxs[v][j] - S.component_frag[c]];
            denom *= M_LN2;
            *grad -= S.weight_matrix->rows[v][j] / denom;
            assert(finite(*grad));
        }

        *grad *= S.tmix[u] + S.tmix[v];
        assert(finite(*grad));
    }

    return p;
}


class MaxPostThread : public InferenceThread
{
    public:
        MaxPostThread(Sampler& S, Queue<ComponentSubset>& q)
            : InferenceThread(S, q)
        {
            topt = nlopt_create(NLOPT_LD_MMA, 1);
            nlopt_set_lower_bounds1(topt, constants::zero_eps);
            nlopt_set_upper_bounds1(topt, 1.0 - constants::zero_eps);
            nlopt_set_xtol_abs1(topt, constants::max_post_x_tolerance);
            nlopt_set_ftol_abs(topt, constants::max_post_objf_tolerance);
            nlopt_set_max_objective(topt, transcript_pair_objf,
                                    reinterpret_cast<void*>(&objfunparam));

            objfunparam.S = &S;
        }

        ~MaxPostThread()
        {
            nlopt_destroy(topt);
        }

        void run_inter_transcript(unsigned int a, unsigned int b);
        void run_inter_component(unsigned int a, unsigned int b);

    private:
        PairwiseObjFuncParam objfunparam;

        /* Optimizer for transcript pair and component pairs, resp. */
        nlopt_opt topt, copt;
};



void MaxPostThread::run_inter_transcript(unsigned int u, unsigned int v)
{
    // XXX: some debug output
    unsigned int c = S.transcript_component[u];
    unsigned int component_size =
        S.component_frag[c + 1] - S.component_frag[c];
    Logger::info("run_inter_transcript: component_size = %u",
                 component_size); 

    if (S.tmix[u] + S.tmix[v] < constants::zero_eps) return;

    objfunparam.u = u;
    objfunparam.v = v;
    double theta = S.tmix[u] / (S.tmix[u] + S.tmix[v]);
    double optf;
    nlopt_result res = nlopt_optimize(topt, &theta, &optf);
    if (res < 0) {
        Logger::warn("Optimization of transcript %u and %u failed.", u, v);
    }
    double z = S.tmix[u] + S.tmix[v];
    S.tmix[u] = theta * z;
    S.tmix[v] = (1.0 - theta) * z;
}


void MaxPostThread::run_inter_component(unsigned int u, unsigned int v)
{
    UNUSED(u);
    UNUSED(v);
    // TODO
}


#if 0
/* Abstract base class for sampler threads. This is specialized in
 * MaxPostThread, and MCMCThread, which estimate maximum posterior and draw
 * samples respectively.
 * */
class SamplerThread
{
    public:
        /* Constructor.
         *
         * Args:
         *   num_componets: Number of connected components.
         *   component_num_transcripts: Number of transcripts in each component.
         *   component_transcripts: An array of transcript indexes for each
         *                          component.
         *   weight_matrix: sparse transcript/fragment weight matrix.
         *   cmix: Component mixture coefficients.
         *   tmix: Within-component transcript mixture coefficients.
         */;
        SamplerThread(Sampler& S,
                      Queue<ComponentSubset>& q)
            : S(S)
            , q(q)
            , thread(NULL)
        {
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            unsigned long seed = reinterpret_cast<unsigned long>(this) *
                                 (unsigned long) time(NULL);
            gsl_rng_set(rng, seed);
        }

        ~SamplerThread()
        {
            gsl_rng_free(rng);
        }

        void run()
        {
            ComponentSubset cs;
            while (true) {
                cs = q.pop();
                if (cs.is_end_of_queue()) break;

                for (unsigned int k = cs.i; k <= cs.j; ++k) {
                    run_intra_component(k);
                }

                for (unsigned int k = cs.i; k <= cs.j - 1; ++k) {
                    run_components(k, k + 1);
                }
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&SamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

        /* Process a single component. */
        void run_intra_component(unsigned int u)
        {
            if (S.component_num_transcripts[u] == 0) return;

            gsl_ran_shuffle(rng, S.component_transcripts[u], 
                            S.component_num_transcripts[u],
                            sizeof(unsigned int));
            for (unsigned int i = 0; i < S.component_num_transcripts[u] - 1; ++i) {
                run_transcripts(S.component_transcripts[u][i],
                                S.component_transcripts[u][i + 1]);
            }
        }

        /* Process a pair of transcripts. */
        virtual void run_transcripts(unsigned int u, unsigned int v) = 0;

        /* Process a pair of a components. */
        virtual void run_components(unsigned int u, unsigned int v) = 0;;

    protected:
        /* To set new mixtures for u and v, tmix'[u] and tmix'[v], we introduce
         * a new variable t, where tmix'[u] = t * (tmix[u] + tmix[v]) and
         * tmix'[v] = (1 - t) * (tmix[u] + tmix[v]).
         *
         * This function computes the derivative of t at the current mix point.
         */
        float transcript_shift_derivative(unsigned int u, unsigned int v)
        {
            unsigned int c = S.transcript_component[u];
            assert(c == S.transcript_component[v]);

            unsigned int component_size =
                S.component_frag[c + 1] - S.component_frag[c];
            float dt = 0.0;
            for (unsigned int j = 0; j < S.weight_matrix->rowlens[u]; ++j) {
                float denom = 
                    S.frag_probs[c][
                        S.weight_matrix->idxs[u][j] - S.component_frag[c]];
                denom *= M_LN2;

                dt += S.weight_matrix->rows[u][j] / denom;
            }

            for (unsigned int j = 0; j < S.weight_matrix->rowlens[v]; ++j) {
                float denom = 
                    S.frag_probs[c][
                        S.weight_matrix->idxs[v][j] - S.component_frag[c]];
                denom *= M_LN2;

                dt -= S.weight_matrix->rows[v][j] / denom;
            }

            dt *= S.tmix[u] + S.tmix[v];

            return finite(dt) ? dt : 0.0;
        }

    private:
        Sampler& S;
        gsl_rng* rng;
        Queue<ComponentSubset>& q;
        boost::thread* thread;
};
#endif


#if 0
/* A thread that works on finding the maximum a posteriori isoform mixture. */
class MaxPostThread
{
    public:
        MaxPostThread(Sampler& S,
                      Queue<unsigned int>& q)
            : S(S)
            , q(q)
            , thread(NULL)
        {
        }

        void run();

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&MaxPostThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

        static const char* task_name;


    private:
        Sampler& S;
        Queue<unsigned int>& q;
        boost::thread* thread;
};


const char* MaxPostThread::task_name = "Optimizing posterior probability";


/* Constrain xs to be on a simplex for nlopt. */
static double simplex_constraint(unsigned int n, const double* xs,
                                 double* grad, void* params)
{
    UNUSED(params);
    double z = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
        z += xs[i];
    }
    if (grad != NULL) std::fill(grad, grad + n, 1.0);

    return z - 1.0;
}


struct ObjFuncParam
{
    unsigned int c; /* component */
    Sampler* S;
};


double transcript_posterior_objf(unsigned int n, const double* xs,
                                 double* grad, void* param)
{
    ObjFuncParam* objfuncparam =
        reinterpret_cast<ObjFuncParam*>(param);
    unsigned int c = objfuncparam->c;
    Sampler& S = *objfuncparam->S;

    unsigned component_size = S.component_frag[c + 1] - S.component_frag[c];
    std::fill(S.frag_probs[c], S.frag_probs[c] + component_size, 0.0f);

    for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
        unsigned int t = S.component_transcripts[c][i];
        asxpy(S.frag_probs[c],
              S.weight_matrix->rows[t],
              xs[i],
              S.weight_matrix->idxs[t],
              S.component_frag[c],
              S.weight_matrix->rowlens[t]);
    }

    double p = dotlog(S.frag_counts[c], S.frag_probs[c], component_size);

    /* extra term to force a unique solution. */
    for (unsigned int i = 0; i < n; ++i) {
        p += xs[i];
    }

    for (unsigned int i = 0; i < component_size; ++i) {
        assert(S.frag_probs[c][i] > 0.0);
        assert(S.frag_probs[c][i] <= 1.0);
        assert(finite(S.frag_probs[c][i]));
    }

    if (grad != NULL) {
        unsigned int off = S.component_frag[c];
        std::fill(grad, grad + n, 0.0);
        for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
            unsigned int t = S.component_transcripts[c][i];
            for (unsigned int j = 0; j < S.weight_matrix->rowlens[t]; ++j) {
                unsigned int k = S.weight_matrix->idxs[t][j];
                grad[i] += S.frag_counts[c][k - off] * S.weight_matrix->rows[t][j] /
                                S.frag_probs[c][k - off];
            }
            grad[i] /= M_LN2;
            assert(grad[i] > 0.0);
            assert(finite(grad[i]));

            /* extra term to force a unique solution. */
            grad[i] += 1.0;
        }

    }

    return p;
}


double component_posterior_objf(unsigned int num_components, const double* cmix,
                                double* grad, void* param)
{
    Sampler& S = *reinterpret_cast<Sampler*>(param);

    double p = 0.0;
    for (unsigned int c = 0; c < num_components; ++c) {
        unsigned component_size = S.component_frag[c + 1] - S.component_frag[c];
        p += dotlogc(S.frag_counts[c], S.frag_probs[c], component_size, cmix[c]);
    }

    if (grad != NULL) {
        std::fill(grad, grad + num_components, 0.0);
        for (unsigned int c = 0; c < num_components; ++c) {
            grad[c] = S.frag_count_sums[c] / (M_LN2 * cmix[c]);
        }
    }

    return p;
}


void MaxPostThread::run()
{
    unsigned int c;
    ObjFuncParam param;
    param.S = &S;
    while (true) {
        c = q.pop();
        if (c >= S.num_components) break;

        unsigned int component_size =
            S.component_frag[c + 1] - S.component_frag[c];

        if (S.component_num_transcripts[c] <= 1 || component_size == 0) {
            /* Set a reasonable starting point. */
            double z = 0.0;
            for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
                S.tmix[S.component_transcripts[c][i]] =
                    (double) S.weight_matrix->rowlens[S.component_transcripts[c][i]] +
                    constants::zero_eps;
                z += S.tmix[S.component_transcripts[c][i]];
            }

            for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
                S.tmix[S.component_transcripts[c][i]] /= z;
            }

            Logger::get_task(task_name).inc();
            continue;
        }

        param.c = c;

        nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP,
                S.component_num_transcripts[c]);
        nlopt_add_inequality_constraint(opt, simplex_constraint,
                                        NULL, constants::simplex_constraint_eps);
        nlopt_set_xtol_abs1(opt, constants::max_post_x_tolerance);
        nlopt_set_ftol_abs(opt, constants::max_post_objf_tolerance);
        nlopt_set_lower_bounds1(opt, constants::zero_eps);
        nlopt_set_upper_bounds1(opt, 1.0);
        nlopt_set_max_objective(opt, transcript_posterior_objf,
                reinterpret_cast<void*>(&param));

        double* xs = new double [S.component_num_transcripts[c]];
        std::fill(xs, xs + S.component_num_transcripts[c],
                1.0 / (double) S.component_num_transcripts[c]);
        double optf;
        nlopt_result res = nlopt_optimize(opt, xs, &optf);
        if (res <= 0) {
            Logger::warn("Posterior probability optimization of component %d failed with code %d.",
                         (int) c, (int) res);
        }

        /* Round tiny numbers to zero and re-normalize. */
        double z = 0.0;
        for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
            if (xs[i] <= 2.0 * constants::zero_eps) xs[i] = 0.0;
            z += xs[i];
        }

        for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
            xs[i] /= z;
        }

        for (unsigned int i = 0; i < S.component_num_transcripts[c]; ++i) {
            S.tmix[S.component_transcripts[c][i]] = xs[i];
        }
        delete [] xs;

        nlopt_destroy(opt);
        Logger::get_task(task_name).inc();
    }
}
#endif



void Sampler::run(unsigned int num_samples)
{
    UNUSED(num_samples); // TODO: delete me

    std::fill(tmix, tmix + ts.size(), 1.0 / (double) ts.size());

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long seed = reinterpret_cast<unsigned long>(this) *
                         (unsigned long) time(NULL);
    gsl_rng_set(rng, seed);

    Queue<ComponentSubset> q;
    std::vector<InferenceThread*> threads(constants::num_threads);
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i] = new MaxPostThread(*this, q);
    }

    unsigned int* cs = new unsigned int [num_components];
    for (unsigned int c = 0; c < num_components; ++c) cs[c] = c;

    // TODO: constant.cpp number of rounds
    for (unsigned int r = 0; r < 10; ++r) {
        gsl_ran_shuffle(rng, cs, num_components, sizeof(unsigned int));
        for (size_t i = 0; i < constants::num_threads; ++i) {
            threads[i]->start();
        }

        // TODO: constants.cpp
        unsigned int c;
        for (c = 0; c + 10 < num_components; c += 10)
        {
            q.push(ComponentSubset(cs + c, 10));
        }

        if (c < num_components) {
            q.push(ComponentSubset(cs + c, num_components - c));
        }

        /* Mark the end. */
        for (size_t i = 0; i < constants::num_threads; ++i) {
            q.push(ComponentSubset());
        }

        for (size_t i = 0; i < constants::num_threads; ++i) {
            threads[i]->join();
        }
    }

    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        fprintf(stderr,
                "%s\t%s\t%e\t%e\t%d\n",
                t->gene_id.get().c_str(),
                t->transcript_id.get().c_str(),
                tmix[t->id],
                transcript_weights[t->id],
                (int) transcript_component[t->id]);
    }
}


void Sampler::init_frag_probs()
{
    for (unsigned int i = 0; i < num_components; ++i) {
        unsigned component_size = component_frag[i + 1] - component_frag[i];
        std::fill(frag_probs[i], frag_probs[i] + component_size, 0.0f);
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
}


