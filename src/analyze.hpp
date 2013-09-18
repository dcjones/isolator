
#ifndef ISOLATOR_ANALYZE_HPP
#define ISOLATOR_ANALYZE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cstdio>

#include "fragment_model.hpp"
#include "queue.hpp"
#include "sampler.hpp"
#include "transcripts.hpp"


class SamplerTickThread;
class TgroupMuSigmaSamplerThread;
class ExperimentTgroupMuSigmaSamplerThread;
class AlphaSampler;
class BetaSampler;


class Analyze
{
    public:
		Analyze(size_t burnin,
                size_t num_samples,
                TranscriptSet& ts,
                const char* genome_filename,
                bool run_gc_correction);
		~Analyze();

		// Add a replicate under a particular condition
		void add_sample(const char* condition_name,
                        const char* filename);

        void run();

    private:
        void setup();
        void cleanup();
        void warmup();
        void sample();

        void qsampler_update_hyperparameters();

        void compute_ts();
        void compute_xs();

        void choose_initial_values();
        void compute_ts_scaling();

        // number of burnin samples
        size_t burnin;

        // number of samples to generate
        size_t num_samples;

        // transcript set
		TranscriptSet& transcripts;

        // File name of a fasta file containing the reference genome sequence
        // against which the reads are aligned.
        const char* genome_filename;

        // True if post-hoc GC content correction should be used.
        bool run_gc_correction;

        // file names for the BAM/SAM file corresponding to each
        std::vector<std::string> filenames;

		// condition index to sample indexes
		std::vector<std::vector<int> > condition_samples;

        // fragment models for each sample
        std::vector<FragmentModel*> fms;

        // quantification samplers for each sample
        std::vector<Sampler*> qsamplers;

        // threads used for iterating samplers
        std::vector<SamplerTickThread*> qsampler_threads;
        std::vector<TgroupMuSigmaSamplerThread*> musigma_sampler_threads;
        std::vector<ExperimentTgroupMuSigmaSamplerThread*> experiment_musigma_sampler_threads;
        AlphaSampler* alpha_sampler;
        BetaSampler* beta_sampler;

        // queues to send work to sampler threads, and be notified on completion
        // of ticks.
        Queue<int> qsampler_tick_queue, qsampler_notify_queue,
                   musigma_sampler_tick_queue, musigma_sampler_notify_queue,
                   experiment_musigma_sampler_tick_queue,
                   experiment_musigma_sampler_notify_queue;


        // matrix containing relative transcript abundance samples, indexed by:
        //   sample -> transcript (tid)
        boost::numeric::ublas::matrix<double> Q;

        // tgroup log-abundance, indexed by sample -> tgroup
        boost::numeric::ublas::matrix<double> ts;

        // transcript abundance relative to other transcripts within the same
        // tgroup. Indexed by replicate -> tid.
        boost::numeric::ublas::matrix<double> xs;

        // tgroup position parameter, indexed by condition -> tgroup
        boost::numeric::ublas::matrix<double> tgroup_mu;

        // tgroup scale parameter, indexed by tgroup
        std::vector<double> tgroup_sigma;

        // parameters of the inverse gamma prior on tgroup_sigma
        double tgroup_alpha, tgroup_beta;

        // experiment-wise tgroup position paremeter, indexed by tgroup
        std::vector<double> experiment_tgroup_mu;

        // experiment-wide tgroup scale parameter, indexed by tgroup
        std::vector<double> experiment_tgroup_sigma;

        // parameters for normal prior over experiment_tgroup_mu
        double experiment_tgroup_mu0, experiment_tgroup_sigma0;

        // experiment wide parameters for the inverse gamma prior on experiment_tgroup_mu
        double experiment_tgroup_alpha, experiment_tgroup_beta;

        // Temporary space used by compute_xs
        std::vector<double> tgroup_expr;

        // Condition index corresponding to the given name
        std::map<std::string, int> condition_index;

        // Condition index of sample i
        std::vector<int> condition;

        // normalization constant for each sample
        std::vector<double> scale;

        // number of sequenced samples
        unsigned int K;

        // number of conditions
        unsigned int C;

        // number of transcripts
        unsigned int N;

        // number of transcription groups (TSS groups, typicall)
        unsigned int T;

        // Hyperparams for inverse gamma prior on tgroup_alpha/tgroup_beta
        double tgroup_nu;
        double tgroup_alpha_alpha, tgroup_beta_alpha;
        double tgroup_alpha_beta, tgroup_beta_beta;
};


#endif

