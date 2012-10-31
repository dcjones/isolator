
#include "constants.hpp"
#include "fragment_model.hpp"
#include "logger.hpp"
#include "queue.hpp"


static const char* param_est_task_name = "Estimating model parameters.";


/* An interval used for fragment model parameter estimation. */
class FragmentModelInterval : public SamScanInterval
{
    public:
        enum IntervalType
        {
            INTERGENIC,
            EXONIC
        } type;

        FragmentModelInterval(const Interval& interval,
                              IntervalType type,
                              Queue<FragmentModelInterval*>& q)
            : SamScanInterval(interval)
            , type(type)
            , q(q)
        {
        }

        /* Once reads have been accumulated in this interval, it adds itself to
         * a work queue to be processed thes deleted. */
        void finish()
        {
            q.push(this);
        }

    private:
        Queue<FragmentModelInterval*>& q;
};


/* A thread used to accumulate statistics to estimate the parameters of
 * the fragment model. */
class FragmentModelThread
{
    public:
        FragmentModelThread(Queue<FragmentModelInterval*>& q)
            : q(q), t(NULL)
        {
        }

        ~FragmentModelThread()
        {
            delete t;
        }

        void run()
        {
            FragmentModelInterval* interval;

            while (true) {
                interval = q.pop();
                if (interval == NULL) break;

                /* dispatch! */
                switch (interval->type) {
                    case FragmentModelInterval::INTERGENIC:
                        /* TODO: measure additive noise */
                        break;

                    case FragmentModelInterval::EXONIC:
                        /* TODO:
                         * measure_fragment_lengths
                         * measure_strand_bias
                         * */

                        /* TODO:
                         * we need to reword the bias correction stuff to be
                         * able to train it here.
                         */
                        break;
                }

                interval->clear();
                Logger::get_task(param_est_task_name).inc();
            }
        }

        void start()
        {
            if (t != NULL) return;
            t = new boost::thread(boost::bind(&FragmentModelThread::run, this));
        }

        void join()
        {
            t->join();
            delete t;
            t = NULL;
        }

    private:
        Queue<FragmentModelInterval*>& q;
        boost::thread* t;
};


FragmentModel::FragmentModel()
{
}


FragmentModel::~FragmentModel()
{
}


void FragmentModel::estimate(TranscriptSet& ts,
                             const char* bam_fn,
                             const char* fa_fn)
{
    if (fa_fn != NULL) {
        /* TODO: train seqbias */
    }

    Queue<FragmentModelInterval*> q(constants::max_estimate_queue_size);

    std::vector<SamScanInterval*> intervals;

    std::vector<Interval> intergenic, exonic;
    ts.get_intergenic(intergenic);
    ts.get_consensus_exonic(exonic);

    /* TODO: if speed is an issue, we can easily subset intergenic and exonic.
     * This is way more data than we really need. */
    std::vector<Interval>::iterator interval;
    for (interval = intergenic.begin(); interval != intergenic.end(); ++interval) {
        intervals.push_back(new FragmentModelInterval(
                            *interval, FragmentModelInterval::INTERGENIC, q));
    }

    for (interval = exonic.begin(); interval != exonic.end(); ++interval) {
        intervals.push_back(new FragmentModelInterval(
                            *interval, FragmentModelInterval::INTERGENIC, q));
    }

    Logger::push_task(param_est_task_name, intervals.size());

    std::vector<FragmentModelThread*> threads(constants::num_threads);
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragmentModelThread(q));
        threads.back()->start();
    }

    sam_scan(intervals, alncnt, bam_fn);

    for (size_t i = 0; i < constants::num_threads; ++i) {
        q.push(NULL);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
        delete threads[i];
    }

    Logger::pop_task(param_est_task_name);
}


