
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

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


typedef std::pair<unsigned int, unsigned int> MultireadFrag;


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
    while (samread(bam_f, b) >= 0) {
        ++read_num;

        if (read_num % 1000 == 0) {
            file_pos = samtell(bam_f);
            if (file_pos >= last_file_pos + input_block_size && input_size > 0) {
                Logger::get_task(task_name).inc();
                last_file_pos = file_pos;
            }
        }

        if (b->core.flag & BAM_FUNMAP ||
            (b->core.flag & BAM_FPAIRED && !(b->core.flag & BAM_FPROPER_PAIR)) ||
            b->core.tid < 0) continue;

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
                TSVec<MultireadFrag>& multiread_frags,
                float* transcript_weights,
                Queue<SamplerInitInterval*>& q)
            : weight_matrix(weight_matrix)
            , frag_counts(frag_counts)
            , multiread_frags(multiread_frags)
            , transcript_weights(transcript_weights)
            , fm(fm)
            , read_indexer(read_indexer)
            , q(q)
            , thread(NULL)
        {
            if (fm.frag_len_dist) {
                frag_len_dist = new EmpDist(*fm.frag_len_dist);
            }
            else {
                frag_len_dist = NULL;
            }
        }

        ~SamplerInitThread()
        {
            delete frag_len_dist;
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
        TSVec<MultireadFrag>& multiread_frags;
        float* transcript_weights;

        FragmentModel& fm;
        Indexer& read_indexer;
        Queue<SamplerInitInterval*>& q;
        boost::thread* thread;

        /* Copy the fragment length distribution to avoid contention between
         * threads. */
        EmpDist* frag_len_dist;

        /* Temprorary space for computing sequence bias, indexed by strand. */
        std::vector<float> mate1_seqbias[2];
        std::vector<float> mate2_seqbias[2];

        /* Exonic length of the transcript whos bias is stored in
         * mate1/mate2_seqbias. */
        pos_t tlen;
        float tw;

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
            // XXX
            //else if (t->transcript_id == "ENST00000233143") {
                //fprintf(stderr, "AAA");
            //}
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
        transcript_sequence_bias(*locus, *t);
        transcript_weights[t->id] = transcript_weight(*t);

        /* If a transcript is very unlikely to have been sequenced at all,
         * ignore it. Otherwise quite a few outliers will be generated in which
         * the annotated transcript is much shorter that what is actually being
         * transcribed. */
        if (transcript_weights[t->id] < constants::min_transcript_weight) {
            continue;
        }

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
        for (j = i; j != k; ++j) {
            sum_weight += j->frag_weight * j->align_pr;
        }
        if (sum_weight == 0.0) continue;

        unsigned int frag_idx = read_indexer.get();
        multiread_frags.push_back(std::make_pair(i->multiread_num, frag_idx));

        /* Sum weight of alignments on the same transcript. */
        float w = 0.0;
        std::vector<MultireadEntry>::iterator j0;
        for (j0 = j = i; j != k; ++j) {
            if (j->transcript_idx != j0->transcript_idx) {
                if (w > 0.0) {
                    weight_matrix.push(j0->transcript_idx, frag_idx, w);
                }
                j0 = j;
                w = 0.0;
            }

            w += j->align_pr * j->frag_weight;
        }
        if (w > 0.0) {
            weight_matrix.push(j0->transcript_idx, frag_idx, w);
        }
    }
}


void SamplerInitThread::transcript_sequence_bias(
                const SamplerInitInterval& locus,
                const Transcript& t)
{
    tlen = t.exonic_length();
    if ((size_t) tlen > mate1_seqbias[0].size()) {
        mate1_seqbias[0].resize(tlen);
        mate1_seqbias[1].resize(tlen);
        mate2_seqbias[0].resize(tlen);
        mate2_seqbias[1].resize(tlen);
    }

    std::fill(mate1_seqbias[0].begin(), mate1_seqbias[0].begin() + tlen, 1.0);
    std::fill(mate1_seqbias[1].begin(), mate1_seqbias[1].begin() + tlen, 1.0);
    std::fill(mate2_seqbias[0].begin(), mate2_seqbias[0].begin() + tlen, 1.0);
    std::fill(mate2_seqbias[1].begin(), mate2_seqbias[1].begin() + tlen, 1.0);

    if (fm.sb == NULL || locus.seq == NULL) return;

    t.get_sequence(tseq0, *locus.seq, fm.sb->getL(), fm.sb->getR());
    t.get_sequence(tseq1, *locus.seq, fm.sb->getR(), fm.sb->getL());
    tseq1.revcomp();

    for (pos_t pos = 0; pos < tlen; ++pos) {
        mate1_seqbias[0][pos] = fm.sb->get_mate1_bias(tseq0, pos + fm.sb->getL());
        mate1_seqbias[1][pos] = fm.sb->get_mate1_bias(tseq1, pos + fm.sb->getL());
        mate2_seqbias[0][pos] = fm.sb->get_mate2_bias(tseq0, pos + fm.sb->getL());
        mate2_seqbias[1][pos] = fm.sb->get_mate2_bias(tseq1, pos + fm.sb->getL());
    }

    std::reverse(mate1_seqbias[1].begin(), mate1_seqbias[1].begin() + tlen);
    std::reverse(mate2_seqbias[1].begin(), mate2_seqbias[1].begin() + tlen);

    if (tlen <= constants::tp_pad) return;

    size_t bin = constants::tp_num_length_bins - 1;
    for (; bin > 0 && tlen < constants::tp_length_bins[bin]; --bin);

    if (t.strand == strand_pos) {
        int d;
        for (pos_t pos = 0; pos < tlen; ++pos) {
            d = constants::tp_num_bins * (pos / (double) (tlen - constants::tp_pad));
            d = std::max<int>(0, std::min<int>(d, constants::tp_num_bins - 1));
            mate1_seqbias[0][pos] *= fm.tp_dist[bin][0][d];

            d = constants::tp_num_bins * ((pos - constants::tp_pad) / (double) (tlen - constants::tp_pad));
            d = std::max<int>(0, std::min<int>(d, constants::tp_num_bins - 1));
            mate1_seqbias[1][pos] *= fm.tp_dist[bin][1][d];
        }
    }
    else {
        int d;
        for (pos_t pos = 0; pos < tlen; ++pos) {
            pos_t p = tlen - pos - 1;
            d = constants::tp_num_bins * ((p - constants::tp_pad) / (double) (tlen - constants::tp_pad));
            d = std::max<int>(0, std::min<int>(d, constants::tp_num_bins - 1));
            mate1_seqbias[0][pos] *= fm.tp_dist[bin][1][d];

            d = constants::tp_num_bins * (p / (double) (tlen - constants::tp_pad));
            d = std::max<int>(0, std::min<int>(d, constants::tp_num_bins - 1));
            mate1_seqbias[1][pos] *= fm.tp_dist[bin][0][d];
        }
    }

    if (t.transcript_id == "ENST00000284292") {

        FILE* out = fopen("ENST00000284292.0.tsv", "w");
        fprintf(out, "i\tw\n");
        for (pos_t i = 0; i < tlen; ++i) {
            fprintf(out, "%ld\t%e\n", i, mate1_seqbias[0][i]);
        }
        fclose(out);

        out = fopen("ENST00000284292.1.tsv", "w");
        fprintf(out, "i\tw\n");
        for (pos_t i = 0; i < tlen; ++i) {
            fprintf(out, "%ld\t%e\n", i, mate1_seqbias[1][i]);
        }
        fclose(out);
    }
}



float SamplerInitThread::transcript_weight(const Transcript& t)
{
    pos_t trans_len = t.exonic_length();
    if ((size_t) trans_len + 1 > ws.size()) ws.resize(trans_len + 1);

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

        /* TODO: The following logic assumes the library type is FR. We need
         * seperate cases to properly handle other library types, that
         * presumably exist somewhere. */

        ws[frag_len] = 0.0;

        float sp0 = t.strand == strand_pos ? fm.strand_specificity :
                                            1.0 - fm.strand_specificity;

        float sp1 = t.strand == strand_neg ? fm.strand_specificity :
                                            1.0 - fm.strand_specificity;

        for (pos_t pos = 0; pos <= trans_len - frag_len; ++pos) {
            // Positive strand fragment weight
            ws[frag_len] +=
                sp0 *
                mate1_seqbias[0][pos] * mate2_seqbias[1][pos + frag_len - 1];

            // Negative strand fragment weight
            ws[frag_len] +=
                sp1 *
                mate2_seqbias[0][pos] * mate1_seqbias[1][pos + frag_len - 1];
        }
    }

    tw = 0.0;

    for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
        float frag_len_pr = frag_len_p(frag_len);
        tw += frag_len_pr * ws[frag_len];
    }

    // XXX
    if (t.transcript_id == "ENST00000233143") {
        FILE* out = fopen("ENST00000233143.tlen.tsv", "w");
        fprintf(out, "frag_len\tp\n");
        for (pos_t frag_len = 1; frag_len <= trans_len; ++frag_len) {
            fprintf(out, "%ld\t%e\n", frag_len, frag_len_p(frag_len));
        }
        fclose(out);
    }

    if (!finite(tw) || tw <= constants::min_transcript_weight) {
        tw = 0.0;
    }

    return tw;
}


float SamplerInitThread::fragment_weight(const Transcript& t,
                                         const AlignmentPair& a)
{
    pos_t frag_len = a.frag_len(t);
    if (frag_len < 0) return 0.0;
    else if (frag_len == 0) {
        pos_t max_frag_len = std::max(a.mate1->end - t.min_start + 1,
                                      t.max_end - a.mate1->start + 1);
        frag_len = std::min(max_frag_len, (pos_t) round(fm.frag_len_med()));
    }

    float w = 1.0;
    if (a.mate1) {
        pos_t offset = t.get_offset(a.mate1->strand == strand_pos ?
                                    a.mate1->start : a.mate1->end);
        if (offset < 0 || offset >= tlen) return 0.0;

        //double cdf;
        //if (a.mate1->strand == t.strand) {
            //cdf = tp_bias[constants::transcript_3p_num_bins - 1][0]->cdf(tlen);
        //}
        //else {
            //cdf = tp_bias[constants::transcript_3p_num_bins - 1][1]->cdf(tlen);
        //}


        w *= mate1_seqbias[a.mate1->strand][offset];
    }

    if (a.mate2) {
        pos_t offset = t.get_offset(a.mate2->strand == strand_pos ?
                                    a.mate2->start : a.mate2->end);
        if (offset < 0 || offset >= tlen) return 0.0;

        w *= mate2_seqbias[a.mate2->strand][offset];
    }

    //pos_t start = t.get_offset(
            //a.mate2 ? std::min<pos_t>(a.mate1->start, a.mate2->start) : a.mate1->start);
    //float c = frag_len_c(tlen - start);

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
    if (frag_len_pr < constants::min_frag_len_pr) return 0.0;
    if (frag_len_pr * w < constants::min_frag_weight) return 0.0;

    //return frag_len_pr / tw;
    return frag_len_pr * w / ws[frag_len] / frag_len_c(tlen);;
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
    Indexer read_indexer;

    /* Collect intervals */
    std::vector<SamplerInitInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    weight_matrix = new WeightMatrix(ts.size());
    transcript_weights = new float [ts.size()];
    TSVec<FragIdxCount> nz_frag_counts;
    TSVec<MultireadFrag> multiread_frags;

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(
                    fm,
                    read_indexer,
                    *weight_matrix,
                    nz_frag_counts,
                    multiread_frags,
                    transcript_weights,
                    q));
        threads.back()->start();
    }

    Logger::debug("Loci: %lu", (unsigned long) intervals.size());

    sam_scan(intervals, bam_fn, fa_fn, "Estimating fragment weights");

    for (size_t i = 0; i < constants::num_threads; ++i) q.push(NULL);
    for (size_t i = 0; i < constants::num_threads; ++i) threads[i]->join();

    /* Free a little space. */
    fm.multireads.clear();

    unsigned int* idxmap = weight_matrix->compact();
    Logger::info("Weight-matrix dimensions: %lu x %lu",
            (unsigned long) weight_matrix->nrow,
            (unsigned long) weight_matrix->ncol);

    /* Update frag_count indexes */
    for (TSVec<FragIdxCount>::iterator i = nz_frag_counts.begin();
            i != nz_frag_counts.end(); ++i) {
        i->first = idxmap[i->first];
    }

    /* Update multiread fragment indexes */
    for (TSVec<MultireadFrag>::iterator i = multiread_frags.begin();
            i != multiread_frags.end(); ++i) {
        i->second = idxmap[i->second];
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

    Logger::debug("Components: %lu", (unsigned long) num_components);

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
        i->second = idxmap[i->second];
    }
    std::sort(multiread_frags.begin(), multiread_frags.end());

    num_multireads = 0;
    unsigned int num_alignments = 0;
    for (unsigned int i = 0; i < multiread_frags.size(); ) {
        unsigned int j;
        for (j = i; j < multiread_frags.size(); ++j) {
            if (multiread_frags[i].first != multiread_frags[j].first) break;
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
            if (multiread_frags[i].first != multiread_frags[j].first) break;
        }

        if (j - i > 1) {
            multiread_alignments[multiread_num] =
                &multiread_alignment_pool[alignment_num];
            multiread_num_alignments[multiread_num] = j - i;

            for (; i < j; ++i) {
                unsigned int c = ds[multiread_frags[i].second];
                multiread_alignment_pool[alignment_num].prob =
                    &frag_probs[c][multiread_frags[i].second - component_frag[c]];
                multiread_alignment_pool[alignment_num].count =
                    &frag_counts[c][multiread_frags[i].second - component_frag[c]];
                ++alignment_num;
            }
            ++multiread_num;
        }

        i = j;
    }

    delete [] ds;
    delete [] idxmap;

    frag_count_sums = new float [num_components];
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
        afree(frag_probs_prop[i]);
        delete [] component_transcripts[i];
    }
    delete [] component_transcripts;
    delete [] component_num_transcripts;
    delete [] multiread_alignment_pool;
    delete [] frag_counts;
    delete [] frag_count_sums;
    delete [] frag_probs;
    delete [] frag_probs_prop;
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
        virtual void run_inter_transcript(unsigned int u, unsigned int v) = 0;

        /* Inference relative expression of two components a and b. */
        virtual void run_component(unsigned int u) = 0;

        void run_intra_component(unsigned int c)
        {
            if (S.component_num_transcripts[c] <= 1) return;

            acopy(S.frag_probs_prop[c], S.frag_probs[c],
                  (S.component_frag[c + 1] - S.component_frag[c]) * sizeof(float));

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

                for (unsigned int i = 0; i < s.len; ++i) {
                    run_component(s.cs[i]);
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


class MCMCThread : public InferenceThread
{
    public:
        MCMCThread(Sampler& S,
                   Queue<ComponentSubset>& q)
            : InferenceThread(S, q)
            , hillclimb(false)
        {
        }

        void run_inter_transcript(unsigned int u, unsigned int v);
        void run_component(unsigned int u);

        /* When true, perform stochastic hill-climbing rather than MCMC. */
        bool hillclimb;

    private:
        float recompute_component_probability(unsigned int u, unsigned int v,
                                              float tmixu, float tmixv,
                                              unsigned int f0, unsigned int f1,
                                              float pf01, float p0, float *d);

        //float recompute_component_probability(unsigned int u, unsigned int v,
                                              //float tmixu, float tmixv,
                                              //unsigned int f0, unsigned int f1,
                                              //float pf01, float p0);

        float transcript_slice_sample_search(float slice_height,
                                             unsigned int u, unsigned int v,
                                             float z0, float p0, bool left);
};


float MCMCThread::recompute_component_probability(unsigned int u, unsigned int v,
                                                  float tmixu, float tmixv,
                                                  unsigned int f0, unsigned int f1,
                                                  float pf01, float p0, float* d)
{
    unsigned int c = S.transcript_component[u];
    assert(c == S.transcript_component[v]);

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

    *d = 0.0;
    *d += asxtydsz(S.frag_counts[c], S.weight_matrix->rows[u], S.frag_probs_prop[c],
                   S.weight_matrix->idxs[u], S.component_frag[c],
                   S.weight_matrix->rowlens[u]);

    *d -= asxtydsz(S.frag_counts[c], S.weight_matrix->rows[v], S.frag_probs_prop[c],
                   S.weight_matrix->idxs[v], S.component_frag[c],
                   S.weight_matrix->rowlens[v]);


    *d *= tmixu + tmixv;
    *d /= M_LN2;

    return p0 - pf01 + dotlog(S.frag_counts[c] + f0, S.frag_probs_prop[c] + f0, f1 - f0);
}


/* An experimental implementation using newton's method. */
float MCMCThread::transcript_slice_sample_search(float slice_height,
                                                 unsigned int u, unsigned int v,
                                                 float z0, float p0, bool left)
{
    static const float peps = 1e-2f;
    static const float zeps = 1e-6f;
    static const float deps = 1e-4f;

    if (fabs(p0 - slice_height) <= peps) return z0;

    /* Find the range of fragment indexes that are affected by updating the u/v
     * mixture, so we can recompute as little as possible. */
    unsigned int comp = S.transcript_component[u];
    unsigned int f0 = S.component_frag[comp+1];
    unsigned int f1 = S.component_frag[comp];
    if (S.weight_matrix->rowlens[u] > 0) {
        f0 = S.weight_matrix->idxs[u][0];
        f1 = 1 + S.weight_matrix->idxs[u][S.weight_matrix->rowlens[u] - 1];
    }

    if (S.weight_matrix->rowlens[v] > 0) {
        f0 = std::min(f0, S.weight_matrix->idxs[v][0]);
        f1 = std::max(f1,  1 + S.weight_matrix->idxs[v][S.weight_matrix->rowlens[v] - 1]);
    }

    if (f0 > f1) f0 = f1 = S.component_frag[comp];
    f0 -= S.component_frag[comp];
    f1 -= S.component_frag[comp];

    /* Round f0 down to the nearst 32-byte boundry so we can use avx/sse. */
    f0 &= 0xfffffff8;

    /* Current probability over [f0, f1). */
    float pf01 = dotlog(S.frag_counts[comp] + f0, S.frag_probs[comp] + f0, f1 - f0);

    const float tmixuv = S.tmix[u] + S.tmix[v];

    /* current estimate */
    float z = left ? constants::zero_eps : 1.0f - constants::zero_eps;

    /* log probability minus slice height */
    float p;

    /* derivative */
    float d;

    z = z0;
    p = recompute_component_probability(
            u, v, z * tmixuv, (1.0f - z) * tmixuv, f0, f1, pf01, p0, &d);
    p -= slice_height;

    /* upper or lower bound on z (depending on 'left') */
    float z_bound = z0;

    unsigned int iter = 0;
    while (fabs(p) > peps && fabs(p/d) > deps) {
        float z1 = z - p / d;

        // Handle the case when the zero is out of bounds.
        if (left && z <= zeps && (z1 < z || p > 0.0)) break;
        if (!left && z >= 1.0 - zeps && (z1 > z || p > 0.0)) break;

        if (p > 0) z_bound = z;

        // Try to void converging to the wrong zero, use bisection when newton
        // tells us to go too far.
        if (left) {
            if (z1 > z0) {
                if (p > 0) {
                    z /= 2.0;
                }
                else {
                    z = (z + z_bound) / 2.0;
                }
            }
            else if (z1 < 0.0) {
                z /= 2.0;
            }
            else z = z1;
        }
        else {
            if (z1 < z0) {
                if (p > 0) {
                    z = (z + 1.0) / 2.0;
                }
                else {
                    z = (z + z_bound) / 2.0;
                }
            }
            else if (z1 > 1.0) {
                z = (z + 1.0) / 2.0;
            }
            else z = z1;
        }

        p = recompute_component_probability(
                u, v, z * tmixuv, (1.0f - z) * tmixuv, f0, f1, pf01, p0, &d);
        p -= slice_height;

        if (++iter > constants::max_newton_iter) break;
    }

    return z;
}


void MCMCThread::run_inter_transcript(unsigned int u, unsigned int v)
{
    unsigned int c = S.transcript_component[u];
    assert(c == S.transcript_component[v]);

    if (S.cmix[c] * (S.tmix[u] + S.tmix[v]) < constants::zero_eps) return;

    unsigned int component_size = S.component_frag[c + 1] - S.component_frag[c];

    double p0 = dotlog(S.frag_counts[c], S.frag_probs[c], component_size);
    float slice_height = hillclimb ? p0 : p0 + fastlog2(gsl_rng_uniform(rng));
    float z0 = S.tmix[u] / (S.tmix[u] + S.tmix[v]);

    float z, s0, s1;
    if (finite(slice_height)) {
        s0 = transcript_slice_sample_search(slice_height, u, v, z0, p0, true);
        s1 = transcript_slice_sample_search(slice_height, u, v, z0, p0, false);
        if (s1 - s0 < constants::zero_eps || s0 > z0 || s1 < z0) {
            return;
        }
        float r = gsl_rng_uniform(rng);
        z = s0 + r * s1 - r * s0;
    }
    else {
        z = gsl_rng_uniform(rng);
    }

    /* proposed tmix[u], tmix[v] values */
    float tmixu = z * (S.tmix[u] + S.tmix[v]);
    float tmixv = (1.0 - z) * (S.tmix[u] + S.tmix[v]);
    float tmix_delta_u = tmixu - S.tmix[u];
    float tmix_delta_v = tmixv - S.tmix[v];

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

    S.tmix[u] = tmixu;
    S.tmix[v] = tmixv;
}


void MCMCThread::run_component(unsigned int u)
{
    float prec = S.frag_count_sums[u] +
                 S.component_num_transcripts[u] * constants::tmix_prior_prec;
    S.cmix[u] = hillclimb ? prec : gsl_ran_gamma(rng, prec, 1.0);
}


/* Two numbers given a range of multireads, representing a chunk of work for the
 * sampler. */
struct MultireadBlock
{
    MultireadBlock()
        : u(INT_MAX)
        , v(INT_MAX)
    {
    }

    MultireadBlock(unsigned int u, unsigned int v)
        : u(u), v(v)
    {
    }

    unsigned int u, v;

    bool is_end_of_queue() const
    {
        return u == INT_MAX;
    }
};


class MultireadSamplerThread
{
    public:
        MultireadSamplerThread(Sampler& S,
                               Queue<MultireadBlock>& q)
            : hillclimb(false)
            , S(S)
            , q(q)
            , thread(NULL)
        {
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            unsigned long seed = reinterpret_cast<unsigned long>(this) *
                                 (unsigned long) time(NULL);
            gsl_rng_set(rng, seed);
        }

        ~MultireadSamplerThread()
        {
            gsl_rng_free(rng);
            if (thread) delete thread;
        }

        void run()
        {
            MultireadBlock block;
            while (true) {
                block = q.pop();
                if (block.is_end_of_queue()) break;

                unsigned int* pc = std::lower_bound(S.component_frag,
                                                    S.component_frag + S.num_components,
                                                    block.u);
                unsigned int c = pc - S.component_frag;

                for (; block.u < block.v; ++block.u) {
                    while (S.component_frag[c + 1] < block.u) ++c;

                    float sumprob = 0.0;
                    unsigned int k = S.multiread_num_alignments[block.u];
                    for (unsigned int i = 0; i < k; ++i) {
                        sumprob += *S.multiread_alignments[block.u][i].prob * S.cmix[c];
                        *S.multiread_alignments[block.u][i].count = 0;
                    }

                    if (hillclimb) {
                        for (unsigned int i = 0; i < k; ++i) {
                            *S.multiread_alignments[block.u][i].count =
                                *S.multiread_alignments[block.u][i].prob * S.cmix[c] / sumprob;
                        }
                    }
                    else {
                        float r = sumprob * gsl_rng_uniform(rng);
                        unsigned int i;
                        for (i = 0; i < k; ++i) {
                            if (r <= *S.multiread_alignments[block.u][i].prob * S.cmix[c]) {
                                break;
                            }
                            else {
                                r -= *S.multiread_alignments[block.u][i].prob * S.cmix[c];
                            }
                        }
                        i = std::min(i, k - 1);
                        *S.multiread_alignments[block.u][i].count = 1;
                    }
                }
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
        boost::thread* thread;
        gsl_rng* rng;
};


void Sampler::run(unsigned int num_samples, SampleDB& out)
{


    /* Initial mixtures */
    // Alternative initialization that prefers the longest trascript.
#if 0
    for (unsigned int i = 0; i < num_components; ++i) {
        unsigned int j_max = 0;
        float max_weight = 0.0;
        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            if (transcript_weights[component_transcripts[i][j]] > max_weight) {
                j_max = j;
                max_weight = transcript_weights[component_transcripts[i][j]];
            }
        }

        tmix[component_transcripts[i][j_max]] =
            1.0 - component_num_transcripts[i] * constants::zero_eps;

        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            if (j != j_max) {
                tmix[component_transcripts[i][j]] = constants::zero_eps;
            }
        }
    }
#endif

#if 0
    for (unsigned int i = 0; i < num_components; ++i) {
        float w = 0.0;
        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            w += transcript_weights[component_transcripts[i][j]];
        }

        for (unsigned int j = 0; j < component_num_transcripts[i]; ++j) {
            tmix[component_transcripts[i][j]] =
                transcript_weights[component_transcripts[i][j]] / w;
        }
    }
#endif

    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        tmix[i] = 1.0 / (float) component_num_transcripts[transcript_component[i]];
    }


    std::fill(cmix, cmix + num_components, 1.0 / (float) num_components);
    init_frag_probs();

    float* maxpost_tmix = new float [weight_matrix->nrow];
    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        maxpost_tmix[i] = tmix[i] * cmix[transcript_component[i]];
    }

    samples = new float * [weight_matrix->nrow];
    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        samples[i] = new float [num_samples];
    }

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long seed = reinterpret_cast<unsigned long>(this) *
                         (unsigned long) time(NULL);
    gsl_rng_set(rng, seed);

    Queue<ComponentSubset> component_queue;
    std::vector<MCMCThread*> mcmc_threads(constants::num_threads);
    for (size_t i = 0; i < constants::num_threads; ++i) {
        mcmc_threads[i] = new MCMCThread(*this, component_queue);
        mcmc_threads[i]->hillclimb = false;
    }

    Queue<MultireadBlock> multiread_queue;
    std::vector<MultireadSamplerThread*> multiread_threads(constants::num_threads);
    for (size_t i = 0; i < constants::num_threads; ++i) {
        multiread_threads[i] = new MultireadSamplerThread(*this, multiread_queue);
        multiread_threads[i]->hillclimb = false;
    }

    unsigned int* cs = new unsigned int [num_components];
    for (unsigned int c = 0; c < num_components; ++c) cs[c] = c;


    unsigned int burnin_samples = constants::sampler_burnin_samples;

    const char* task_name = "Sampling";
    Logger::push_task(task_name, num_samples + burnin_samples +
                                 constants::sampler_hillclimb_samples);

    for (unsigned int sample_num = 0;
            sample_num < num_samples + constants::sampler_hillclimb_samples; ) {
        /* Sample multiread alignments */
        for (size_t i = 0; i < constants::num_threads; ++i) {
            multiread_threads[i]->start();
        }

        unsigned int r;
        for (r = 0; r + constants::sampler_multiread_block_size < num_multireads;
                r += constants::sampler_multiread_block_size) {
            multiread_queue.push(MultireadBlock(
                        r, r + constants::sampler_multiread_block_size));
        }
        if (r < num_multireads) {
            multiread_queue.push(MultireadBlock(r, num_multireads));
        }

        for (size_t i = 0; i < constants::num_threads; ++i) {
            multiread_queue.push(MultireadBlock());
        }

        for (size_t i = 0; i < constants::num_threads; ++i) {
            multiread_threads[i]->join();
        }

        /* Recompute total fragments in each component */
        std::fill(frag_count_sums, frag_count_sums + num_components, 0.0f);
        for (size_t i = 0; i < num_components; ++i) {
            unsigned int component_size = component_frag[i + 1] - component_frag[i];
            for (unsigned int j = 0; j < component_size; ++j) {
                frag_count_sums[i] += frag_counts[i][j];
            }
        }

        /* Sample transcript abundance. */
        gsl_ran_shuffle(rng, cs, num_components, sizeof(unsigned int));
        for (size_t i = 0; i < constants::num_threads; ++i) {
            mcmc_threads[i]->start();
        }

        unsigned int c;
        for (c = 0; c + constants::sampler_component_block_size < num_components;
                c += constants::sampler_component_block_size)
        {
            component_queue.push(ComponentSubset(
                        cs + c, constants::sampler_component_block_size));
        }
        if (c < num_components) {
            component_queue.push(ComponentSubset(cs + c, num_components - c));
        }

        /* Mark the end. */
        for (size_t i = 0; i < constants::num_threads; ++i) {
            component_queue.push(ComponentSubset());
        }

        for (size_t i = 0; i < constants::num_threads; ++i) {
            mcmc_threads[i]->join();
        }

        /* normalize cmix */
        float z = 0.0;
        for (unsigned int c = 0; c < num_components; ++c) {
            z += cmix[c];
        }

        for (unsigned int c = 0; c < num_components; ++c) {
            cmix[c] /= z;
        }


        // DEBUGING XXX
        //for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
            //if (t->transcript_id == "ENST00000233143") {
                //Logger::info("ENST00000233143 has %f fragments.",
                             //frag_count_sums[transcript_component[t->id]]);
            //}
        //}


        /* Record a new sample */
        if(burnin_samples == 0 && sample_num < num_samples) {
            double total_weight = 0.0, new_total_weight = 0.0;
            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                if (transcript_weights[i] < constants::min_transcript_weight) {
                    samples[i][sample_num] = 0.0;
                }
                else {
                    samples[i][sample_num] =
                        tmix[i] * cmix[transcript_component[i]] /
                        transcript_weights[i];

                    total_weight += samples[i][sample_num];
                }
            }

            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                samples[i][sample_num] /= total_weight;

                if (samples[i][sample_num] < constants::round_down_eps) {
                    samples[i][sample_num] = 0.0;
                }
                new_total_weight += samples[i][sample_num];
            }

            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                samples[i][sample_num] /= new_total_weight;
            }

            ++sample_num;

            if (sample_num == num_samples) {
                for (size_t i = 0; i < constants::num_threads; ++i) {
                    mcmc_threads[i]->hillclimb = true;
                    //multiread_threads[i]->hillclimb = true;
                }
            }
        }
        /* Update the current maximum posterior */
        else if (burnin_samples == 0) {
            double total_weight = 0.0, new_total_weight = 0.0;
            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                if (transcript_weights[i] < constants::min_transcript_weight) {
                    maxpost_tmix[i] = 0.0;
                }
                else {
                    maxpost_tmix[i] =
                        tmix[i] * cmix[transcript_component[i]] /
                        transcript_weights[i];

                    total_weight += maxpost_tmix[i];
                }
            }

            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                maxpost_tmix[i] /= total_weight;

                if (maxpost_tmix[i] < constants::round_down_eps) {
                    maxpost_tmix[i] = 0.0;
                }
                new_total_weight += maxpost_tmix[i];
            }

            for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
                maxpost_tmix[i] /= new_total_weight;
            }

            ++sample_num;
        }
        /* Throw a sample away */
        else {
            --burnin_samples;
        }

        /* Check for numerical errors. */
        for (unsigned int i = 0; i < num_components; ++i) {
            for (unsigned int j = 0; j < component_frag[i + 1] - component_frag[i]; ++j ) {
                if (!finite(frag_probs[i][j])) {
                    Logger::warn("numerical error: %f", frag_probs[i][j]);
                }
            }
        }

        Logger::get_task(task_name).inc();
    }

    Logger::pop_task(task_name);

    out.begin_transaction();
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        out.insert_sampler_result(
                t->transcript_id, t->gene_id,
                transcript_weights[t->id],
                maxpost_tmix[t->id],
                samples[t->id],
                num_samples);
    }
    out.commit_transaction();

    for (unsigned int i = 0; i < weight_matrix->nrow; ++i) {
        delete [] samples[i];
    }
    delete [] samples;
    delete [] maxpost_tmix;
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


