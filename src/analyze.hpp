
#ifndef ISOLATOR_ANALYZE_HPP
#define ISOLATOR_ANALYZE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cstdio>

#include "hdf5.hpp"
#include "fragment_model.hpp"
#include "queue.hpp"
#include "sampler.hpp"
#include "transcripts.hpp"


class SamplerTickThread;
class ConditionTgroupMuSigmaSamplerThread;
class ExperimentTgroupMuSigmaSamplerThread;
class GammaBetaSampler;
class AlphaSampler;
class BetaSampler;
class GammaNormalSigmaSampler;
class ConditionSpliceMuSigmaEtaSamplerThread;
class ExperimentSpliceMuSigmaSamplerThread;
typedef std::pair<int, int> IdxRange;


class Analyze
{
    public:
        Analyze(unsigned int rng_seed,
                size_t burnin,
                size_t num_samples,
                TranscriptSet& ts,
                const char* genome_filename,
                bool run_gc_correction,
                bool run_3p_correction,
                double experiment_tgroup_sigma_alpha,
                double experiment_tgroup_sigma_beta,
                double experiment_splice_sigma_alpha,
                double experiment_splice_sigma_beta,
                double condition_tgroup_alpha,
                double condition_tgroup_beta_a,
                double condition_tgroup_beta_b,
                double condition_splice_alpha,
                double condition_splice_beta_a,
                double condition_splice_beta_b);
        ~Analyze();

        // Add a replicate under a particular condition
        void add_sample(const char* condition_name,
                        const char* filename);

        void run(hid_t file_id, bool dryrun);

    private:
        void setup_samplers();
        void setup_output(hid_t output_file_id);
        void cleanup();
        void warmup();
        void sample(bool optimize_state);
        void write_output(size_t sample_num);

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

        // True if GC content correction should be used.
        bool run_gc_correction;

        // True if 3' bias should be corrected
        bool run_3p_correction;

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
        std::vector<ConditionTgroupMuSigmaSamplerThread*> musigma_sampler_threads;
        std::vector<ExperimentTgroupMuSigmaSamplerThread*> experiment_musigma_sampler_threads;
        GammaBetaSampler* gamma_beta_sampler;
        BetaSampler* invgamma_beta_sampler;
        GammaNormalSigmaSampler* gamma_normal_sigma_sampler;
        std::vector<ConditionSpliceMuSigmaEtaSamplerThread*> splice_mu_sigma_sampler_threads;
        std::vector<ExperimentSpliceMuSigmaSamplerThread*>
            experiment_splice_mu_sigma_sampler_threads;

        // queues to send work to sampler threads, and be notified on completion
        // of ticks.
        Queue<int> qsampler_tick_queue, qsampler_notify_queue;

        // work is doled out in block for these. Otherwise threads can starve
        // when there are few sample in the experiment
        Queue<IdxRange> musigma_sampler_tick_queue,
                        experiment_musigma_sampler_tick_queue,
                        splice_mu_sigma_sampler_tick_queue,
                        experiment_splice_mu_sigma_sampler_tick_queue;

        Queue<int> musigma_sampler_notify_queue,
                   experiment_musigma_sampler_notify_queue,
                   splice_mu_sigma_sampler_notify_queue,
                   experiment_splice_mu_sigma_sampler_notify_queue;

        // We maintain a different rng for every unit of work for threads.
        // That way we can actually make isolator run reproducible.
        std::vector<rng_t> splice_rng_pool;
        std::vector<rng_t> tgroup_rng_pool;

        // matrix containing relative transcript abundance samples, indexed by:
        //   sample -> transcript (tid)
        boost::numeric::ublas::matrix<double> Q;

        // tgroup log-abundance, indexed by sample -> tgroup
        boost::numeric::ublas::matrix<double> ts;

        // transcript abundance relative to other transcripts within the same
        // tgroup. Indexed by replicate -> tid.
        boost::numeric::ublas::matrix<double> xs;

        // tgroup position parameter, indexed by condition -> tgroup
        boost::numeric::ublas::matrix<double> condition_tgroup_mu;

        // tgroup scale parameter, indexed by tgroup
        std::vector<double> condition_tgroup_sigma;

        // parameters of the inverse gamma prior on condition_tgroup_sigma
        double condition_tgroup_alpha, condition_tgroup_beta;

        // parameters of the inverse gamma prior on condition_splice_sigma
        double condition_splice_alpha, condition_splice_beta;

        // experiment-wise tgroup position paremeter, indexed by tgroup
        std::vector<double> experiment_tgroup_mu;

        // experiment-wide tgroup scale parameter
        double experiment_tgroup_sigma;

        // gamma hypeparameters for prior on experiment_tgroup_sigma
        double experiment_tgroup_sigma_alpha;
        double experiment_tgroup_sigma_beta;

        // degrees of freedom parameter for experiment tgroup mu prior
        double experiment_tgroup_nu;

        // parameters for normal prior over experiment_tgroup_mu
        double experiment_tgroup_mu0, experiment_tgroup_sigma0;

        // temporary space used by compute_xs
        std::vector<double> tgroup_expr;

        // tids belonging to each tgroup (indexed by tgroup)
        std::vector<std::vector<unsigned int> > tgroup_tids;

        // sorted indexes of tgroups with multiple transcripts
        std::vector<unsigned int> spliced_tgroup_indexes;

        // condition slice mean indexed by condition, spliced tgroup, transcript
        // according to spliced_tgroup_indexes and tgroup_tids
        std::vector<std::vector<std::vector<double> > > condition_splice_mu;

        // per-spliced-tgroup experiment wide logistic-normal mean
        std::vector<std::vector<double> > experiment_splice_mu;

        // prior parameters for experimient_splice_mu
        double experiment_splice_nu, experiment_splice_mu0, experiment_splice_sigma0;

        // experiment std. dev.
        double experiment_splice_sigma;

        // gamma hypeparameters for prior on experiment_splice_sigma
        double experiment_splice_sigma_alpha;
        double experiment_splice_sigma_beta;

        // splicing precision, indexed by spliced tgroup
        std::vector<std::vector<double> > condition_splice_sigma;

        // overparameterization to unstick stuck samplers
        std::vector<std::vector<double> > condition_splice_eta;

        // flattened condition_splice_sigma used for sampling alpha, beta
        // params.
        std::vector<double> condition_splice_sigma_work;
        std::vector<double> experiment_splice_sigma_work;
        std::vector<double> experiment_tgroup_sigma_work;

        // paramaters for the inverse gamma priors on splice_alpha and
        // splice_beta
        double condition_splice_beta_a,
               condition_splice_beta_b;

        // Condition index corresponding to the given name
        std::map<std::string, int> condition_index;

        // Condition index of sample i
        std::vector<int> condition;

        // normalization constant for each sample
        std::vector<double> scale;

        // temporary space used for computing scale
        std::vector<double> scale_work;

        // number of sequenced samples
        unsigned int K;

        // number of conditions
        unsigned int C;

        // number of transcripts
        unsigned int N;

        // number of transcription groups (TSS groups, typicall)
        unsigned int T;

        // Hyperparams for inverse gamma prior on tgroup_alpha/tgroup_beta
        double condition_tgroup_beta_a,
               condition_tgroup_beta_b;

        // RNG used for alpha/beta samplers
        unsigned int rng_seed;
        rng_t rng;

        // HDF5 dataspace ids, for output purposes

        // dataspaces
        hid_t h5_experiment_tgroup_dataspace;
        hid_t h5_condition_tgroup_dataspace;
        hid_t h5_tgroup_row_mem_dataspace;
        hid_t h5_sample_quant_dataspace;
        hid_t h5_sample_quant_mem_dataspace;
        hid_t h5_experiment_splice_dataspace;
        hid_t h5_condition_splice_mu_dataspace;
        hid_t h5_condition_splice_sigma_dataspace;
        hid_t h5_splicing_mem_dataspace;
        hid_t h5_sample_scaling_dataspace;
        hid_t h5_sample_scaling_mem_dataspace;

        // datasets
        hid_t h5_experiment_mean_dataset;
        hid_t h5_experiment_sd_dataset;
        hid_t h5_condition_mean_dataset;
        hid_t h5_sample_quant_dataset;
        hid_t h5_experiment_splice_mu_dataset;
        hid_t h5_experiment_splice_sigma_dataset;
        hid_t h5_condition_splice_mu_dataset;
        hid_t h5_condition_splice_sigma_dataset;
        hid_t h5_sample_scaling_dataset;

        // variable length array for splicing paramaters
        hid_t h5_splice_param_type;

        // structure for the ragged-array splicing data
        hvl_t* h5_splice_work;

        // a write buffer for converting doubles to floats before hdf5 output
        std::vector<float> tgroup_row_data;

        friend void write_qc_data(FILE* fout, Analyze& analyze);
};


#endif

