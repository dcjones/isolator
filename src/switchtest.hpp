
#ifndef ISOLATOR_SWITCHTEST_HPP
#define ISOLATOR_SWITCHTEST_HPP

#include <map>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_rng.h>
#include "intervals.hpp"
#include "sample_db.hpp"
#include "transcripts.hpp"


class LambdaSliceSampler;


class SwitchTest
{
	public:
		SwitchTest(TranscriptSet& ts);
		~SwitchTest();

		// Add a replicate under a particular condition
		void add_replicate(const char* condition_name, SampleDB& replicate_data);

		// Run the sampler.
		void run(unsigned int num_samples);

	private:
		void check_data();
		void build_sample_matrix();

		// compute tss usage for replicate i, sample k and store it in the appropriate
		// row in 'tss_usage'
		void compute_tss_usage(unsigned int i, unsigned int k);

		// compute the log-likelihood of sample k within replicate i, given mu and sigma.
		double sample_log_likelihood(unsigned int i, unsigned int k);

		void sample_replicate_transcript_abundance();
		void sample_condition_tss_usage();

		// number of replicates
		unsigned int n;

		// number of transcripts
		unsigned int m;

		typedef std::map<Interval, std::vector<unsigned int> > TSMap;
		typedef std::pair<const Interval, std::vector<unsigned int> > TSMapItem;
		TSMap tss_group;
		TSMap tss_tts_group;

		// map each transcript index to it's tss index
		std::vector<unsigned int> tss_index;

		TranscriptSet& ts;

		// map a condition name to data from some number of replicates
		typedef std::map<std::string, std::vector<SampleDB*> > CondMap;
		typedef std::pair<const std::string, std::vector<SampleDB*> > CondMapItem;
		CondMap data;

		// condition index to condition name
		std::vector<std::string> condition_name;

		// condition index to replicate indexes
		std::vector<std::vector<unsigned int> > condition_replicates;

		// replicate index to condition index
		std::vector<unsigned int> replicate_condition;

		// temporary space used to marginalize tss abundance, indexed by
		// replicate index -> tss index
		boost::numeric::ublas::matrix<float> tss_usage;

		// matrices containing sampled transcript abundances from 'isolator quantify'
		std::vector<boost::numeric::ublas::matrix<float> > samples;

		// replicate transcript abundance (index into samples)
		std::vector<unsigned int> repl_sample_idx;

		// condition tss abundance parameters mean and precision
		boost::numeric::ublas::matrix<float> mu;
		boost::numeric::ublas::matrix<float> lambda;

		// parameters for the normal prior on mu, indexd by
		// replicate index -> tss index
		double mu0, lambda0;

		// gamma parameters for the prior on sigma
		std::vector<float> alpha, beta;

		LambdaSliceSampler* lambda_sampler;

		gsl_rng* rng;
};


#endif
