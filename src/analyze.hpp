
#ifndef ISOLATOR_ANALYZE_HPP
#define ISOLATOR_ANALYZE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "fragment_model.hpp"
#include "queue.hpp"
#include "sampler.hpp"
#include "transcripts.hpp"


class SamplerTickThread;

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

        void qsampler_tick();
        void qsampler_update_hyperparameters(const std::vector<double>& params);

        // XXX: Fuuuuck! I have two things named ts: tgroup expression and
        // "transcript set".
        //
        // Options:
        //   Rename the transcript set to "transcripts". I should really do this
        //   everywhere if I'm going to do it.

        void compute_ts(std::vector<std::vector<double> >& ts);
        void compute_xs(std::vector<std::vector<double> >& xs);

        void choose_initial_values(std::vector<double>& cont_params,
                                   std::vector<int>& disc_params);

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
		std::vector<std::vector<unsigned int> > condition_samples;

        // fragment models for each sample
        std::vector<FragmentModel*> fms;

        // quantification samplers for each sample
        std::vector<Sampler*> qsamplers;

        // threads used for iterating samplers
        std::vector<SamplerTickThread*> qsampler_threads;

        // queues to send work to sampler threads, and be notified on completion
        // of ticks.
        Queue<int> qsampler_tick_queue;
        Queue<int> qsampler_tock_queue;

        // matrix containing relative transcript abundance samples, indexed by:
        //   replicate -> transcript (tid)
        boost::numeric::ublas::matrix<double> Q;

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

        // used to load data into the sampler
        friend class AnalyzeSamplerData;
};


#endif

