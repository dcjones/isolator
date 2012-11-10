
#include "constants.hpp"
#include "logger.hpp"
#include "queue.hpp"
#include "sam_scan.hpp"
#include "sampler.hpp"

class SamplerInitInterval : public SamScanInterval
{
    public:
        SamplerInitInterval(
                TranscriptSetLocus locus,
                Queue<SamplerInitInterval*>& q)
            : locus(locus)
            , q(q)
        {
        }

        void finish()
        {
            q.push(this);
        }


    private:
        TranscriptSetLocus locus;
        Queue<SamplerInitInterval*>& q;

        friend class SamplerInitThread;
};


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
            SamplerInitInterval* interval;
            while (true) {
                if ((interval = q.pop()) == NULL) break;

                /* Collect loci. */

                /* TODO:
                 * The interval contains a set of transcripts along
                 * with the reads than inersect that locus.
                 *
                 * We need to produce from that a matrix of fragment
                 * likelihoods.
                 */
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
        FragmentModel& fm;
        Queue<SamplerInitInterval*>& q;
        boost::thread* thread;
};



Sampler::Sampler(const char* bam_fn, const char* ref_fn,
                 TranscriptSet& ts,
                 FragmentModel& fm)
    : ts(ts)
    , fm(fm)
{
    /* Producer/consumer queue of intervals containing indexed reads to be
     * processed. */
    Queue<SamplerInitInterval*> q;

    /* Collect intervals */
    std::vector<SamScanInterval*> intervals;
    for (TranscriptSetLocusIterator i(ts); i != TranscriptSetLocusIterator(); ++i) {
        intervals.push_back(new SamplerInitInterval(*i, q));
    }

    std::vector<SamplerInitThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new SamplerInitThread(fm, q));
        threads.back()->start();
    }

    /* TODO: Run sam_scan, or some version of.
     *
     * We can't use the regular sam_scan function, we let it get too
     * specialized.
     *
     * Possibilities:
     *   1. rework sam_scan
     *   2. just write a seperate version.
     *
     *
     * If we choose 1, it seems like the current version should just become part
     * of fragment_model.cpp...
     */

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
