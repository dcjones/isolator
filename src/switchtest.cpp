
#include <set>
#include <gsl/gsl_randist.h>
#include "constants.hpp"
#include "switchtest.hpp"
#include "logger.hpp"


using namespace boost::numeric::ublas;


static double sq(double x)
{
	return x * x;
}

static const double LN_SQRT_TWO_PI = log(sqrt(2.0 * M_PI));

static double gaussian_lnpdf(double mu, double lambda, double x)
{
	return -LN_SQRT_TWO_PI + sqrt(lambda) - sq(x - mu) * lambda  / 2.0;
}


SwitchTest::SwitchTest(TranscriptSet& ts)
	: ts(ts)
	, mu0(-10.0)
	, lambda0(0.1)
{
	// organize transcripts by TSS and TTS
	for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
		pos_t tss_pos = t->strand == strand_pos ? t->min_start : t->max_end;
		Interval tss(t->seqname, tss_pos, tss_pos, t->strand);
		tss_group[tss].push_back(t->id);

		Interval tss_tts(t->seqname, t->min_start, t->max_end, t->strand);
		tss_tts_group[tss_tts].push_back(t->id);
	}

	m = ts.size();
	tss_index.resize(m);
	unsigned int tss_id = 0;
	for (std::map<Interval, std::vector<unsigned int> >::iterator i = tss_group.begin();
		 i != tss_group.end(); ++i, ++tss_id) {
		for (std::vector<unsigned int>::iterator j = i->second.begin();
			 j != i->second.end(); ++j) {
			tss_index[*j] = tss_id;
		}
	}

    rng = gsl_rng_alloc(gsl_rng_mt19937);
}


SwitchTest::~SwitchTest()
{
	gsl_rng_free(rng);
}


void SwitchTest::add_replicate(const char* condition_name, SampleDB& replicate_data)
{
	data[condition_name].push_back(&replicate_data);
}


void SwitchTest::run(unsigned int num_samples)
{
	check_data();
	build_sample_matrix();

	n = 0;
	for (CondMap::iterator i = data.begin(); i != data.end(); ++i) {
		n += i->second.size();
	}

	mu.resize(n, tss_group.size());
	lambda.resize(n, tss_tts_group.size());
	tss_usage.resize(n, tss_group.size());
	// TODO: initialize

	repl_sample_idx.resize(n);
	std::fill(repl_sample_idx.begin(), repl_sample_idx.end(), 0);


	// TODO burnin

	for (unsigned int sample_num = 0; sample_num < num_samples; ++sample_num) {
		sample_replicate_transcript_abundance();
		sample_condition_tss_usage();
	}

	// TODO:
	// We need to add parameter vectors for
	//   replicate tss usage
	//   replicate tts usage
	//   replatete splice usage

	//   condition tss usage
	//   condition tts usage
	//   condition splice usage

	//   hyperparameters for priors on
	//     variability of


	samples.clear();
}


// basic sanity check for input data: it must have the same transcript_ids and be non-empty.
void SwitchTest::check_data()
{
	if (data.size() == 0) {
		Logger::abort("No data to test.");
	}
	else if (data.size() == 1) {
		Logger::abort("Only one condition provided to test.");
	}

	std::set<TranscriptID> transcript_ids;
	std::set<TranscriptID> replicate_transcript_ids;
	unsigned int repl_num = 0;
	for (CondMap::iterator i = data.begin(); i != data.end(); ++i) {
		if (i->second.size() == 0) {
			Logger::abort("Condition %s has no replicates.", i->first.c_str());
		}

		replicate_transcript_ids.clear();
		std::vector<SampleDB*>::iterator j;
		for (j = i->second.begin(); j != i->second.end(); ++j, ++repl_num) {
			for (SampleDBIterator k(**j); k != SampleDBIterator(); ++k) {
				replicate_transcript_ids.insert(k->transcript_id);
			}
		}

		if (repl_num == 0) {
			transcript_ids = replicate_transcript_ids;
		}
		else if (transcript_ids != replicate_transcript_ids) {
			Logger::abort("Same was generated using different gene annotations. "
				          "Use the same annotations to proceed.");
		}
	}
}


void SwitchTest::build_sample_matrix()
{
	condition_name.clear();
	condition_replicates.clear();

	std::map<TranscriptID, unsigned int> tids;
	for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
		tids[t->transcript_id] = t->id;
	}

	// for each replicate build a sample matrix in which rows are indexed by transcript
	// and columns by sample number.
	unsigned int repl_num = 0;
	for (CondMap::iterator i = data.begin(); i != data.end(); ++i) {
		condition_name.push_back(i->first);
		condition_replicates.push_back(std::vector<unsigned int>());
		std::vector<SampleDB*>::iterator j;
		for (j = i->second.begin(); j != i->second.end(); ++j, ++repl_num) {
			condition_replicates.back().push_back(repl_num);
			samples.push_back(matrix<float>());
			samples.back().resize(tids.size(), (*j)->get_num_samples());
			for (SampleDBIterator k(**j); k != SampleDBIterator(); ++k) {
				for (unsigned int l = 0; l < k->samples.size(); ++l) {
					samples.back()(tids[k->transcript_id], l) = k->samples[l];
				}
			}
		}
	}
}


void SwitchTest::compute_tss_usage(unsigned int i, unsigned int k)
{
	std::fill(tss_usage.data().begin(), tss_usage.data().end(), 0.0);
	for (unsigned int j = 0; j < m; ++j) {
		tss_usage(i, tss_index[j]) += samples[i](j, k);
	}
}


double SwitchTest::sample_log_likelihood(unsigned int i, unsigned int k)
{
	double p = 0.0;

	// tss usage log-probabilities
	compute_tss_usage(i, k);
	for (unsigned int j = 0; j < tss_usage.size2(); ++j)  {
		double x = std::max<double>(tss_usage(i, j), constants::zero_eps);
		p += gaussian_lnpdf(mu(i, j), lambda(i, j), x);
	}

	// tts usage log-probabilities
	// TODO:
	// For every replicate and tss, we need a vector of probabilities for
	// the catagorial distribution. Then, another level for tss/tts/splicing.
	// Fancy shit.

	// splicing log-probabilities
	// TODO

	return p;
}


void SwitchTest::sample_replicate_transcript_abundance()
{
	for (unsigned int i = 0; i < samples.size(); ++i) {
		// current sample log probability
		double q = sample_log_likelihood(i, repl_sample_idx[i]);

		unsigned int proposal = gsl_rng_uniform_int(rng, samples[i].size2());
		double p = sample_log_likelihood(i, proposal);

		// metropolis accept/reject
		if (log(gsl_rng_uniform(rng)) < p - q) {
			repl_sample_idx[i] = proposal;
		}
		else {
			compute_tss_usage(i, repl_sample_idx[i]);
		}
	}
}


void SwitchTest::sample_condition_tss_usage()
{
	for (unsigned int i = 0; i < condition_replicates.size(); ++i) {
		for (unsigned int k = 0; k < tss_usage.size2(); ++k) {

			// TODO: Draw precision sample
			// We have to sample precision numerically
			// Should we fall back on our old standard of slice sampling?


			// Ok, this is what we need to tackle tomorrow morning.
			// This is a comparatively easy slice sampler. 

#if 0
			unsigned num_cond_replicates = condition_replicates[i].size();
			double a = alpha + (double) num_cond_replicates / 2.0;

			double sample_mu = 0.0;
			for (std::vector<unsigned int>::iterator j = condition_replicates[i].begin();
				 j != condition_replicates[i].end(); ++j) {
				sample_mu += tss_usage(*j, k);
			}
			sample_mu /= (double) num_cond_replicates;

			double sample_var = 0.0
			for (std::vector<unsigned int>::iterator j = condition_replicates[i].begin();
				 j != condition_replicates[i].end(); ++j) {
				sample_var += sq(sample_mu - tss_usage(*j, k));
			}

			double b = beta +
			           sample_var / 2.0 +
			           num_cond_replicates * sq(mu0 - sample_mean)
			}
#endif


			// sample from mu given lambda and transcript abundance
			unsigned num_cond_replicates = condition_replicates[i].size();
			double sample_mu = 0.0;
			for (std::vector<unsigned int>::iterator j = condition_replicates[i].begin();
				 j != condition_replicates[i].end(); ++j) {
				sample_mu += tss_usage(*j, k);
			}
			sample_mu /= (double) num_cond_replicates;

			double weighted_lambda_ik = num_cond_replicates * lambda(i, k);
			double numer = sample_mu * weighted_lambda_ik + mu0 * lambda0;
			double denom = lambda0 + weighted_lambda_ik;
			double posterior_mu = numer / denom;
			double posterior_sigma = sqrt(1.0 / denom);
			mu(i, k) = posterior_mu + gsl_ran_gaussian(rng, posterior_sigma);
		}
	}
}
