
#ifndef ISOLATOR_ANALYZE_HPP
#define ISOLATOR_ANALYZE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <set>
#include <string>
#include <vector>

#include "fragment_model.hpp"
#include "sampler.hpp"
#include "transcripts.hpp"


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
        //void compute_depth();
        //void choose_initial_values(std::vector<double>& cont_params,
                                   //std::vector<int>& disc_params);

        // number of burnin samples
        size_t burnin;

        // number of samples to generate
        size_t num_samples;

        // transcript set
		TranscriptSet& ts;

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

        // matrix containing relative transcript abundance samples, indexed by:
        //   replicate -> transcript (tid)
        boost::numeric::ublas::matrix<float> Q;

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

        // used to load data into the sampler
        friend class AnalyzeSamplerData;
};


#endif

