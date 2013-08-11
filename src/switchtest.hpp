
#ifndef ISOLATOR_SWITCHTEST_HPP
#define ISOLATOR_SWITCHTEST_HPP

#include <map>
#include <string>
#include <vector>
#include <set>
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_rng.h>
#include "intervals.hpp"
#include "sample_db.hpp"
#include "transcripts.hpp"


class LambdaSliceSampler;
class AlphaSliceSampler;
class BetaSliceSampler;


class SwitchTest
{
	public:
		SwitchTest(TranscriptSet& ts);
		~SwitchTest();

		// Add a replicate under a particular condition
		void add_replicate(const char* condition_name, SampleDB& replicate_data);

		// Run the sampler.
		void run(unsigned int num_samples, const std::vector<double>& quantiles);

	private:
		void check_data();
		void build_sample_matrix();

        // initialize the splicing_alpha and splicing beta arrays
        void init_splicing_lambda_priors();

        // compute tss usage for replicate i, sample k and store it in the
        // appropriate row in 'tss_usage'
		void compute_tss_usage_splicing(unsigned int i, unsigned int k);

		// compute the log-likelihood of sample k within replicate i, given mu and sigma.
		double sample_log_likelihood(unsigned int i, unsigned int k);

		void sample_replicate_transcript_abundance();
		void sample_condition_tss_usage();
        void sample_condition_splicing();

        void output_mu_samples(FILE* out, const std::vector<double>& quantiles,
                               unsigned int cond1, unsigned int cond2);

		// number of replicates
		unsigned int m;

		// number of transcription start sites
		unsigned int n_tss;

		typedef std::map<Interval, std::vector<unsigned int> > TSMap;
		typedef std::pair<const Interval, std::vector<unsigned int> > TSMapItem;
		TSMap tss_group;
		TSMap tss_tts_group;

		// map each transcript index to it's tss index
		std::vector<unsigned int> tss_index;

        // Transcript and gene IDs indexd by tss index.
        std::vector<std::set<GeneID> > tss_gene_ids;
        std::vector<std::set<TranscriptID> > tss_transcript_ids;
        std::vector<std::vector<unsigned int> > tss_tids;

		TranscriptSet& ts;

        // tid to transcript_ids and gene_ids
        std::vector<TranscriptID> transcript_ids;
        std::vector<GeneID> gene_ids;

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

        // base precision for tss_usage accounted for by poisson sampling
		// indexed by: replicate index -> tss index
        boost::numeric::ublas::matrix<float> tss_base_precision;

        // factor by which a replicate's tss_expression should be normalized
        std::vector<double> tss_usage_norm_factor;
        std::vector<double> tss_usage_norm_factor_work;

		// matrices containing sampled transcript abundances from 'isolator quantify'
		// Indexed by: replicate -> transcript -> sample_num
		std::vector<boost::numeric::ublas::matrix<float> > samples;

        // effective lengths, indexed by: replicate -> tid
        boost::numeric::ublas::matrix<float> effective_lengths;

        // numbers of reads, conditioned by replicate
        std::vector<float> read_counts;

		// replicate transcript abundance (index into samples)
		std::vector<unsigned int> repl_sample_idx;

		// condition tss abundance mean. Indexed by: condition -> tss.
		boost::numeric::ublas::matrix<float> mu;

		// tss abundance precision. Indexed by: tss.
		std::vector<float> lambda;

		// gamma parameters for the gamma prior on lambda[i]
		double alpha, beta;

		// gamma parameters for the gamma pior on the alpha
		double alpha_alpha_0, beta_alpha_0;

		// gamma parameters for the gamma pior on the beta
		double alpha_beta_0, beta_beta_0;

		// parameters for the normal prior on mu, indexd by
		// replicate index -> tss index
		double mu0, lambda0;

        // tts and splicing, given tss usage.
        // indexed by replicate index -> transcript id
        boost::numeric::ublas::matrix<float> splicing;

        // condition splicing means, indexd by condition
        boost::numeric::ublas::matrix<float> splicing_mu;

        // tss-wise splicing precisions, indexd by tss index
        std::vector<float> splicing_lambda;

        // pooled gamma prior on splicing precision, conditioned on the number of
        // isoforms in the group.
        // Indexed by: simplex degree (isoforms per gene)
        std::vector<float> splicing_alpha;
        std::vector<float> splicing_beta;

        // symetric dirichlet prior on splicing means
        double splicing_mu0;

        // splicing "hyper-priors"
        double splicing_alpha_alpha_0, splicing_beta_alpha_0;
        double splicing_alpha_beta_0, splicing_beta_beta_0;

		LambdaSliceSampler* lambda_sampler;
		AlphaSliceSampler* alpha_sampler;
		BetaSliceSampler* beta_sampler;

        // collected samples indexed by:
        // sample_num -> condition -> tss
        std::vector<boost::numeric::ublas::matrix<float> > mu_samples;

		gsl_rng* rng;
};


#endif
