
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


/* Store information needed to process multireads. */
struct Multiread
{
    Multiread()
    {
    }

    /* Information necessary to compute one entry in the likelihood matrix. */
    struct Cell
    {
        unsigned int transcript_id;
        float align_pr;
        float frag_pr;
    };

    std::vector<Cell> cells;
};


/* When initializing the likelihood matrix rows for reads with multiple
 * alignments cannot be added until all alignments have been read.
 * This class collects information for multireads to be processed later. */
class MultireadIndex
{
    public:
        MultireadIndex()
        {
            t = hattrie_create();
        }

        ~MultireadIndex()
        {
            hattrie_free(t);
        }

        void push(const char* id, const Multiread::Cell& cell)
        {
            boost::lock_guard<boost::mutex> lock(mut);
            Multiread** r = reinterpret_cast<Multiread**>(
                    hattrie_get(t, id, strlen(id)));

            if (*r == NULL) {
                *r = new Multiread();
            }

            (*r)->cells.push_back(cell);
        }

    private:
        hattrie_t* t;
        boost::mutex mut;
};


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

    private:
        int32_t tid;
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
    faidx_t* fa_f = fai_load(fa_fn);
    if (fa_f == NULL) {
        Logger::abort("Can't open FASTA file %s.", fa_fn);
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
    fai_destroy(fa_f);

    if (task_name) Logger::pop_task(task_name);
}


class SamplerInitThread
{
    public:
        SamplerInitThread(
                FragmentModel& fm,
                Queue<SamplerInitInterval*>& q)
            : fm(fm)
            , q(q)
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

        /* Compute (a number proportional to) the probability of observing the
         * given fragment from the given transcript, assuming every transcript
         * is equally expressed. */
        float fragment_probability(
                const Transcript& t,
                const AlignmentPair& a);

        FragmentModel& fm;
        Queue<SamplerInitInterval*>& q;
        boost::thread* thread;
};


void SamplerInitThread::process_locus(SamplerInitInterval* locus)
{
    /* Number of transcripts */
    size_t n = locus->ts.size();

    /* A map of fragments onto sequential indices. */
    typedef std::map<AlignmentPair, unsigned int> FragIdx;
    FragIdx fragidx;

    /* A set of fragments. */
    typedef std::set<AlignmentPair> FragSet;

    /* The set of fragments within this locus that have no compatible
     * transcripts. */
    FragSet excluded_frags;

    /* Assign indexes to reads, filter that are not compatible with any
     * transcript, and handle multireads. */
    for (ReadSetIterator r(locus->rs); r != ReadSetIterator(); ++r) {
        bool multiread = fm.multireads.has(r->first);
        /* TODO: do something with multireads. */
        if (multiread) continue;

        for (AlignedReadIterator a(*r->second); a != AlignedReadIterator(); ++a) {
            TranscriptSetLocus::iterator t;
            bool has_compatible_transcript;
            for (t = locus->ts.begin(); t != locus->ts.end(); ++t) {
                // TODO: compute fragment likelihood
                // set flag if non-zero
            }

            if (has_compatible_transcript) {
                unsigned int idx = fragidx.size();
                fragidx[*a] = idx;
            }
            else {
                excluded_frags.insert(*a);
            }
        }
    }



    /* TODO */

    /* Collect loci. */

    /* TODO:
     * The interval contains a set of transcripts along
     * with the reads than inersect that locus.
     *
     * We need to produce from that a matrix of fragment
     * likelihoods.
     */


    /* Something like this:
     *
     * If this is not a multiread:
     *      Make a new row for it in the loci's matrix.
     * Otherwise:
     *      ???
     *
     *      Otherwise what?
     *
     *      We have to figure this out before we write a bunch of
     *      code that will end up in the garbage.
     *
     *
     */
}


float SamplerInitThread::fragment_probability(const Transcript& t,
                                              const AlignmentPair& a)
{
    pos_t fragment_length;

    /* The fragment is not compatible with the trascript. */
    if (fragment_length < 0) return 0.0;

    /* The fragment is compatibly, but lacks a mate, so we guess the fragment
     * length. */
    if (fragment_length == 0) {
        // TODO
    }

    return 0.0;
}



Sampler::Sampler(const char* bam_fn, const char* fa_fn,
                 TranscriptSet& ts,
                 FragmentModel& fm)
    : ts(ts)
    , fm(fm)
{
    /* Producer/consumer queue of intervals containing indexed reads to be
     * processed. */
    Queue<SamplerInitInterval*> q;

    /* Collect intervals */
    std::vector<SamplerInitInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(fm, q));
        threads.back()->start();
    }

    sam_scan(intervals, bam_fn, fa_fn, "Initializing samplers");

    for (size_t i = 0; i < constants::num_threads; ++i) {
        q.push(NULL);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
        delete threads[i];
    }

    /* Now, we should have a bunch of rectangular sub-matrices. One for each
     * loci. We need to combine these into a larger matrix, then add rows for
     * multireads somehow. Yikes.
     */
}


Sampler::~Sampler()
{
}
