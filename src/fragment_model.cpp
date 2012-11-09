
#include <boost/unordered_map.hpp>

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


/* A thread used to accumulate statistics to estimate the parameters of
 * the fragment model. */
class FragmentModelThread
{
    public:
        FragmentModelThread(Queue<FragmentModelInterval*>& q)
            : q(q), t(NULL)
        {
            strand_bias[0] = strand_bias[1] = 0;
        }

        ~FragmentModelThread()
        {
            delete t;
        }

        void run()
        {
            FragmentModelInterval* interval;
            ReadSet::UniqueReadCounts read_counts;

            while (true) {
                interval = q.pop();
                if (interval == NULL) break;

                /* dispatch! */
                if (interval->type == FragmentModelInterval::INTERGENIC) {
                    /* TODO: measure additive noise */
                }
                else if (interval->type == FragmentModelInterval::EXONIC) {
                    interval->rs.make_unique_read_counts(read_counts);
                    measure_fragment_lengths(read_counts);
                    measure_strand_bias(interval->strand, read_counts);
                    read_counts.clear();
                }

                interval->clear();
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

        /* strand_bias[0] counts the number of times the read's strand agrees
         * with the transcripts, and strand_bias[1] the number of disagreements.
         * */
        unsigned long strand_bias[2];

        /* A hash table mapping fragment lengths to number of observations */
        boost::unordered_map<unsigned int, unsigned int> frag_lens;

        /* A hash table mapping a number k to the number of intergenic
         * positions with k read starts. */
        boost::unordered_map<unsigned int, unsigned int> noise_counts;

    private:
        void measure_fragment_lengths(const ReadSet::UniqueReadCounts& counts)
        {
            boost::unordered_map<unsigned int, unsigned int>::iterator c;
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (j->mate1 == NULL || j->mate2 == NULL) continue;
                    pos_t len = j->naive_frag_len();
                    if (len <= 0 || len > constants::max_frag_len) continue;

                    c = frag_lens.find(len);
                    if (c == frag_lens.end()) {
                        frag_lens.insert(std::make_pair((unsigned int) len, 1));
                    }
                    else {
                        c->second += 1;
                    }
                }
            }
        }

        void measure_strand_bias(strand_t strand,
                                 const ReadSet::UniqueReadCounts& counts)
        {
            //Logger::info("A");
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                //Logger::info("B");
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (!j->mate1) continue;
                    //Logger::info("C");
                    strand_bias[j->mate1->strand == strand ? 0 : 1]++;
                }
            }
        }

        Queue<FragmentModelInterval*>& q;
        boost::thread* t;
};


FragmentModel::FragmentModel()
    : sb(NULL)
    , frag_len_dist(NULL)
{
}


FragmentModel::~FragmentModel()
{
    if (sb) delete sb;
    if (frag_len_dist) delete frag_len_dist;
}


void FragmentModel::estimate(TranscriptSet& ts,
                             const char* bam_fn,
                             const char* fa_fn)
{
    Queue<FragmentModelInterval*> q(constants::max_estimate_queue_size);

    std::vector<SamScanInterval*> intervals;

    std::vector<Interval> exonic;
    ts.get_consensus_exonic(exonic);
    std::vector<Interval>::iterator interval;
    for (interval = exonic.begin(); interval != exonic.end(); ++interval) {
        intervals.push_back(new FragmentModelInterval(
                                    *interval,
                                    FragmentModelInterval::EXONIC,
                                    q));
    }

    std::vector<FragmentModelThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragmentModelThread(q));
        threads.back()->start();
    }

    PosTable mate1_pos_tab, mate2_pos_tab;
    sam_scan(intervals, alncnt,
             mate1_pos_tab, mate2_pos_tab,
             bam_fn, "Indexing reads");

    if (fa_fn != NULL) {
        sb = new sequencing_bias(fa_fn,
                                 mate1_pos_tab, mate2_pos_tab,
                                 constants::seqbias_num_reads,
                                 constants::seqbias_left_pos,
                                 constants::seqbias_right_pos);
    }
    else sb = NULL;

    for (size_t i = 0; i < constants::num_threads; ++i) {
        q.push(NULL);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
    }

    /* Collect strand specificity statistics. */
    unsigned long strand_bias[2] = {0, 0};
    for (size_t i = 0; i < constants::num_threads; ++i) {
        strand_bias[0] += threads[i]->strand_bias[0];
        strand_bias[1] += threads[i]->strand_bias[1];
    }

    strand_specificity = (float) strand_bias[0];
    strand_specificity /= (float) (strand_bias[0] + strand_bias[1]);

    Logger::info("Strand specificity: %0.1f%%", 100.0 * strand_specificity);

    /* Collect fragment length statistics. */
    std::map<unsigned int, unsigned int> frag_lens;
    boost::unordered_map<unsigned int, unsigned int>::iterator f;
    std::map<unsigned int, unsigned int>::iterator g;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        for (f = threads[i]->frag_lens.begin();
             f != threads[i]->frag_lens.end();
             ++f) {
            g = frag_lens.find(f->first);
            if (g == frag_lens.end()) {
                frag_lens.insert(*f);
            }
            else {
                g->second += f->second;
            }
        }
    }

    if (frag_lens.size() > constants::frag_len_min_pe_reads) {
        unsigned int* frag_len_vals = new unsigned int [frag_lens.size()];
        unsigned int* frag_len_lens = new unsigned int [frag_lens.size()];
        unsigned int sum_lens = 0;

        size_t i;
        for (g = frag_lens.begin(), i = 0; g != frag_lens.end(); ++g) {
            frag_len_vals[i] = g->first;
            frag_len_lens[i] = g->second;
            sum_lens += g->second;
        }

        /* prune fragment lengths with very low probability */
        for (i = 0; i < frag_lens.size(); ++i) {
            if ((float) frag_len_lens[i] / (float) sum_lens <
                    constants::min_frag_len_pr) {
                frag_len_lens[i] = 0;
            }
        }

        frag_len_dist = new EmpDist(frag_len_vals, frag_len_lens, frag_lens.size(),
                                    constants::frag_len_dist_smoothing);

        delete frag_len_vals;
        delete frag_len_lens;

        Logger::info("Fragment length distribution estimated. (Median: %d)",
                     (int) frag_len_dist->median());
    }
    else {
        Logger::info("Too few paired-end reads to estimate fragment length distribution.");
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete threads[i];
    }
}

