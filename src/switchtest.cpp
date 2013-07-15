
#include <set>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
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
	return -LN_SQRT_TWO_PI + sqrt(lambda) - sq(x - mu) * lambda / 2.0;
}


static double gamma_lnpdf(double alpha, double beta, double x)
{
	return alpha * log(beta) -
	       gsl_sf_lngamma(alpha) +
	       (alpha - 1) * log(x) -
	       beta * x;
}


// A generic slice one dimensional sampler.
class SliceSampler
{
public:
	SliceSampler();
	virtual ~SliceSampler();

	// Generate a new sample given the current value x0.
	double sample(double x0);

	// Log-probability at value x.
	virtual double lpr(double x) = 0;

private:
	double find_slice_edge(double x0, double slice_height, double init_step);

	gsl_rng* rng;
};


SliceSampler::SliceSampler()
{
    rng = gsl_rng_alloc(gsl_rng_mt19937);
}


SliceSampler::~SliceSampler()
{
	gsl_rng_free(rng);
}


double SliceSampler::sample(double x0)
{
	double y0 = lpr(x0);
	double slice_height = log(gsl_rng_uniform(rng)) + y0;

	// find slice extents
	double x_min = find_slice_edge(x0, slice_height, 1e-1);
	double x_max = find_slice_edge(x0, slice_height, 1e+1);

	// sample
	double x = 0.0;
	double y = -INFINITY;
	while (y < slice_height) {
		x = x_min + (x_max - x_min) * gsl_rng_uniform(rng);
		y = lpr(y);
	}
	return x;
};


double SliceSampler::find_slice_edge(double x0, double slice_height, double step)
{
	const double eps = 1e-2;

	// step further and further until something outside the slice is found
	double x, y;
	do {
		x = x0 + step;
		y = lpr(x);
		step *= 2;
	} while (y > slice_height);
	step /= 2;

	// binary search to find the edge
	double a = x0, b = x0 + step;
	while (fabs(b - a) < eps) {
		double c = (a + b) / 2;
		double w = lpr(c);

		// No! this depends on the direction
		if (w > slice_height) a = c;
		else                  b = c;
	}

	return (a + b) / 2;
}


// Slice sampler realizations
class LambdaSliceSampler : public SliceSampler
{
public:
	~LambdaSliceSampler() {}

	double lpr(double lambda)
	{
		double p = gamma_lnpdf(alpha, beta, lambda);
		for (std::vector<double>::iterator x = xs.begin(); x != xs.end(); ++x) {
			p += gaussian_lnpdf(mu, lambda, *x);
		}
		return p;
	}

	double alpha;
	double beta;
	double mu;
	std::vector<double> xs;
};



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

	Logger::info("%u transcript start sites", (unsigned int) tss_group.size());

	m = ts.size();
	tss_index.resize(m);
	unsigned int tss_id = 0;
	BOOST_FOREACH (const TSMapItem& i, tss_group) {
		BOOST_FOREACH (const unsigned int& j, i.second) {
			tss_index[j] = tss_id;
		}
		++tss_id;
	}

	lambda_sampler = new LambdaSliceSampler();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
}


SwitchTest::~SwitchTest()
{
	delete lambda_sampler;
	gsl_rng_free(rng);
}


void SwitchTest::add_replicate(const char* condition_name, SampleDB& replicate_data)
{
	data[condition_name].push_back(&replicate_data);
}


void SwitchTest::run(unsigned int num_samples)
{
	n = 0;
	BOOST_FOREACH (CondMapItem& i, data) n += i.second.size();

	build_sample_matrix();

	mu.resize(n, tss_group.size());
	std::fill(mu.data().begin(), mu.data().end(), mu0);

	lambda.resize(n, tss_group.size());
	std::fill(lambda.data().begin(), lambda.data().end(), lambda0);

	alpha.resize(tss_group.size());
	std::fill(alpha.begin(), alpha.end(), 1.0);

	beta.resize(tss_group.size());
	std::fill(beta.begin(), beta.end(), 0.5);

	tss_usage.resize(n, tss_group.size());
	repl_sample_idx.resize(n);
	std::fill(repl_sample_idx.begin(), repl_sample_idx.end(), 0);


	// TODO burnin
	for (unsigned int sample_num = 0; sample_num < num_samples; ++sample_num) {
		sample_replicate_transcript_abundance();
		sample_condition_tss_usage();
	}

	samples.clear();
}


void SwitchTest::build_sample_matrix()
{
	const char* task_name = "Loading quantification data";
	Logger::push_task(task_name, n);

	if (data.size() == 0) {
		Logger::abort("No data to test.");
	}
	else if (data.size() == 1) {
		Logger::abort("Only one condition provided to test.");
	}

	condition_name.clear();
	condition_replicates.clear();

	std::map<TranscriptID, unsigned int> tids;
	for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
		tids[t->transcript_id] = t->id;
	}

	// for each replicate build a sample matrix in which rows are indexed by transcript
	// and columns by sample number.
	unsigned int repl_idx = 0;
	BOOST_FOREACH (CondMapItem& i, data) {
		if (i.second.size() == 0) {
			Logger::abort("Condition %s has no replicates.", i.first.c_str());
		}

		condition_name.push_back(i.first);
		condition_replicates.push_back(std::vector<unsigned int>());
		BOOST_FOREACH (SampleDB*& j, i.second) {
			condition_replicates.back().push_back(repl_idx);
			samples.push_back(matrix<float>());
			samples.back().resize(tids.size(), j->get_num_samples());
			for (SampleDBIterator k(*j); k != SampleDBIterator(); ++k) {
				std::map<TranscriptID, unsigned int>::iterator tid;
				tid = tids.find(k->transcript_id);
				if (tid == tids.end()) {
					Logger::abort(
						"Transcript %s is in the quantification data but "
						" not in the provided annotations.",
						k->transcript_id.get().c_str());
				}

				matrix_row<matrix<float> > row(samples.back(), tid->second);
				std::copy(k->samples.begin(), k->samples.end(), row.begin());
			}
			Logger::get_task(task_name).inc();
			++repl_idx;
		}
	}

	Logger::pop_task(task_name);
}


void SwitchTest::compute_tss_usage(unsigned int i, unsigned int k)
{
	matrix_row<matrix<float> > row(tss_usage, i);
	std::fill(row.begin(), row.end(), 0.0);
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
		p += gaussian_lnpdf(mu(i, j), lambda(i, j), log(x));
	}

	// TODO: tts usage and splicing probabilities
	// For every replicate and tss, we need a vector of probabilities for
	// the catagorial distribution. Then, another level for tss/tts/splicing.
	// Fancy shit.

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
			// sample from lambda given mu and transcript abundance
			lambda_sampler->alpha = alpha[k];
			lambda_sampler->beta = beta[k];
			lambda_sampler->mu = mu(i, k);
			lambda_sampler->xs.resize(condition_replicates[i].size());
			matrix_column<matrix<float> > col(tss_usage, k);
			std::copy(col.begin(), col.end(), lambda_sampler->xs.begin());
			lambda(i, k) = lambda_sampler->sample(lambda(i, k));

			// sample from mu given lambda and transcript abundance
			unsigned num_cond_replicates = condition_replicates[i].size();
			double sample_mu = 0.0;
			BOOST_FOREACH(float& x, col) {
				sample_mu += x;
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
