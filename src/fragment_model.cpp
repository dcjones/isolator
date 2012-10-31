
#include "constants.hpp"
#include "fragment_model.hpp"
#include "logger.hpp"
#include "queue.hpp"


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
    const char* task_name = "Estimating model parameters.";

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

    Logger::push_task(task_name, intervals.size());

    /* TODO: start up a bunch of threads to read from q */

    sam_scan(intervals, alncnt, bam_fn, fa_fn);

    /* TODO: use the collected data to estimate some things. */

    Logger::pop_task(task_name);
}


