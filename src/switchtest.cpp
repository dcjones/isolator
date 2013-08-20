
#include <numeric>
#include <set>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "constants.hpp"
#include "dirichlet.h"
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
    return -LN_SQRT_TWO_PI + log(sqrt(lambda)) - sq(x - mu) * lambda / 2.0;
}


static double gamma_lnpdf(double alpha, double beta, double x)
{
    return alpha * log(beta) -
           gsl_sf_lngamma(alpha) +
           (alpha - 1) * log(x) -
           beta * x;
}


static double ran_scaled_inv_chisq(gsl_rng* rng, double v, double var)
{
    return (var * v) / gsl_ran_chisq(rng, v);
}


// Log-pdf for student's t-distribution, avoiding some expensive recomputation.
class StudentsTLogPdf
{
    public:
        // v: degrees of freedom
        StudentsTLogPdf(double v)
            : v(v)
        {
            base = gsl_sf_lngamma((v + 1)/2)
                 - gsl_sf_lngamma(v/2)
                 - log(sqrt(M_PI * v));
        }

        double operator()(double mu, double lambda, double x) const
        {
            return base
                 - log(mu)
                 - ((v + 1) / 2) * log1p(sq(x - mu) * lambda / v);
        }

    private:
        double base;
        double v;
};


// A generic slice one dimensional sampler.
class SliceSampler
{
public:
    SliceSampler(double lower_limit, double upper_limit);
    virtual ~SliceSampler();


protected:
    // Log-probability at value x.
    virtual double lpr(double x) = 0;

    // Generate a new sample given the current value x0.
    double sample(double x0);

    // Bounds on the parameter being sampled over
    double lower_limit, upper_limit;

    gsl_rng* rng;
private:
    double find_slice_edge(double x0, double slice_height, double init_step);

};


SliceSampler::SliceSampler(double lower_limit, double upper_limit)
    : lower_limit(lower_limit)
    , upper_limit(upper_limit)
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
    double step = 0.1 * x0;
    double x_min = find_slice_edge(x0, slice_height, -step);
    double x_max = find_slice_edge(x0, slice_height, +step);

    // sample
    double x = 0.0;
    double y;
    x = x_min + (x_max - x_min) * gsl_rng_uniform(rng);
    y = lpr(x);

    return x;
};


double SliceSampler::find_slice_edge(double x0, double slice_height, double step)
{
    const double eps = 1e-2 * x0;

    // step further and further until something outside the slice is found
    double x, y;
    do {
        x = std::max(lower_limit, std::min(upper_limit, x0 + step));
        y = lpr(x);
        step *= 2;
    } while (y > slice_height && x > lower_limit && x < upper_limit);
    step /= 2;

    // binary search to find the edge
    double a, b;
    if (step < 0.0) {
        a = std::max(lower_limit, std::min(upper_limit, x0 + step));
        b = x0;
    }
    else {
        a = x0;
        b = std::max(lower_limit, std::min(upper_limit, x0 + step));
    }

    while (fabs(b - a) > eps) {
        double c = (a + b) / 2;
        double w = lpr(c);

        if (step < 0) {
            if (w > slice_height) b = c;
            else                  a = c;
        }
        else {
            if (w > slice_height) a = c;
            else                  b = c;
        }
    }

    return (a + b) / 2;
}


// Slice sampler realizations
class LambdaSliceSampler : public SliceSampler
{
public:
    LambdaSliceSampler(double mu_dof,
                       double& alpha,
                       double& beta,
                       boost::numeric::ublas::matrix<float>& mu,
                       boost::numeric::ublas::matrix<float>& tss_usage,
                       boost::numeric::ublas::matrix<float>& tss_base_precision,
                       std::vector<unsigned int>& replicate_condition)
        : SliceSampler(1e-4, INFINITY)
        , tpdf(mu_dof)
        , alpha(alpha)
        , beta(beta)
        , mu(mu)
        , tss_usage(tss_usage)
        , tss_base_precision(tss_base_precision)
        , replicate_condition(replicate_condition)
    {}

    virtual ~LambdaSliceSampler() {}


    double sample(unsigned int tss_idx, double current_lambda)
    {
        this->tss_idx = tss_idx;
        return SliceSampler::sample(current_lambda);
    }

protected:

    double lpr(double lambda)
    {
        // prior
        double p = gamma_lnpdf(alpha, beta, lambda);

        // likelihood
        for (unsigned int i = 0; i < tss_usage.size1(); ++i) {
            p += gaussian_lnpdf(
                    mu(replicate_condition[i], tss_idx),
                    (lambda * tss_base_precision(i, tss_idx)) /
                    (lambda + tss_base_precision(i, tss_idx)),
                    tss_usage(i, tss_idx));
        }

        return p;
    }

private:
    unsigned tss_idx;

    StudentsTLogPdf tpdf;

    double& alpha;
    double& beta;
    boost::numeric::ublas::matrix<float>& mu;
    boost::numeric::ublas::matrix<float>& tss_usage;
    boost::numeric::ublas::matrix<float>& tss_base_precision;
    std::vector<unsigned int>& replicate_condition;
};


class AlphaSliceSampler : public SliceSampler
{
public:
    AlphaSliceSampler(
            double alpha_alpha_0,
            double beta_alpha_0,
            std::vector<float>& lambda)
        : SliceSampler(1e-3, INFINITY)
        , alpha_alpha_0(alpha_alpha_0)
        , beta_alpha_0(beta_alpha_0)
        , lambda(lambda)
    {}

    virtual ~AlphaSliceSampler() {}

    double sample(double beta, double current_alpha)
    {
        this->beta = beta;
        return SliceSampler::sample(current_alpha);
    }

protected:
    double lpr(double alpha)
    {
        // prior
        double p = gamma_lnpdf(alpha_alpha_0, beta_alpha_0, alpha);

        // likelihood
        BOOST_FOREACH (float& l, lambda) {
            p += gamma_lnpdf(alpha, beta, l);
        }

        return p;
    }

private:
    double beta;
    double alpha_alpha_0, beta_alpha_0;
    std::vector<float>& lambda;
};


class BetaSliceSampler : public SliceSampler
{
public:
    BetaSliceSampler(
            double alpha_beta_0,
            double beta_beta_0,
            std::vector<float>& lambda)
        : SliceSampler(1e-8, INFINITY)
        , alpha_beta_0(alpha_beta_0)
        , beta_beta_0(beta_beta_0)
        , lambda(lambda)
    {}

    virtual ~BetaSliceSampler() {}

    double sample(double alpha, double current_beta)
    {
        this->alpha = alpha;
        return SliceSampler::sample(current_beta);
    }

protected:
    double lpr(double beta)
    {
        // prior
        double p = gamma_lnpdf(alpha_beta_0, beta_beta_0, beta);

        // likelihood
        BOOST_FOREACH (float& l, lambda) {
            p += gamma_lnpdf(alpha, beta, l);
        }

        return p;
    }

private:
    double alpha;
    double alpha_beta_0, beta_beta_0;
    std::vector<float>& lambda;
};


class SplicingMuSliceSampler : public SliceSampler
{
    public:
        SplicingMuSliceSampler(double splicing_mu0,
                               matrix<float>& splicing,
                               matrix<float>& splicing_mu,
                               std::vector<float>& splicing_lambda,
                               std::vector<std::vector<unsigned int> >& tss_tids,
                               std::vector<std::vector<unsigned int> >&
                                    condition_replicates)
            : SliceSampler(0.0, 1.0)
            , splicing_mu0(splicing_mu0)
            , splicing(splicing)
            , splicing_mu(splicing_mu)
            , splicing_lambda(splicing_lambda)
            , tss_tids(tss_tids)
            , condition_replicates(condition_replicates)
        {}

        void sample(unsigned int cond, unsigned int tss_id)
        {
            if (tss_tids[tss_id].size() < 2) return;

            this->tss_id = tss_id;

            for (unsigned int round = 0; round < tss_tids[tss_id].size(); ++round) {
                u = gsl_rng_uniform_int(rng, tss_tids[tss_id].size());
                v = gsl_rng_uniform_int(rng, tss_tids[tss_id].size() - 1);
                if (v >= u) ++v;

                double mu_u = splicing_mu(cond, tss_tids[tss_id][u]),
                       mu_v = splicing_mu(cond, tss_tids[tss_id][v]);

                double x0 = mu_u / (mu_u + mu_v);

                double x = SliceSampler::sample(x0);
                splicing_mu(cond, tss_tids[tss_id][u]) = x * (mu_u + mu_v);
                splicing_mu(cond, tss_tids[tss_id][v]) = (1.0 - x) * (mu_u + mu_v);
            }
        }

    protected:
        // we sample the relative expression of two transcritps at a time. The
        // indexes are held here.
        unsigned int u, v;

        // tss id being sampled over
        unsigned int tss_id;

        // condition being samples over
        unsigned int cond;


        double lpr(double x)
        {
            double prec = splicing_lambda[tss_id];

            // TODO: No! this is wrong. I'm evaluating at x0, not x.[o

            // prior
            double p = 0.0;
            BOOST_FOREACH (unsigned int tid, tss_tids[tss_id]) {
                p += (splicing_mu0 - 1) * log(splicing_mu(cond, tid));
            }
            p -= tss_tids[tss_id].size() * gsl_sf_lngamma(splicing_mu0) -
                 gsl_sf_lngamma(splicing_mu0 * tss_tids[tss_id].size());

            // likelihood
            return p;
        }

    private:
        double splicing_mu0;

        matrix<float>& splicing;
        matrix<float>& splicing_mu;
        std::vector<float>& splicing_lambda;
        std::vector<std::vector<unsigned int> >& tss_tids;
        std::vector<std::vector<unsigned int> >& condition_replicates;
};

// Blah!  Let's concentrate on mu first.
#if 0
class SplicingLambdaSliceSampler : public SliceSampler
{
    public:
        SplicingLambdaSliceSampler()
            : SliceSampler(1e-4, INFINITY)
        {
        }

        virtual ~SplicingLambdaSliceSampler() {}


        double sample()
        {
        }

    protected:

        double lpr(double lambda)
        {
            // prior
            double p = 

        };


    private:
};
#endif

SwitchTest::SwitchTest(TranscriptSet& ts)
    : ts(ts)
    , mu_dof(2.0)
    , alpha_alpha_0(5)
    , beta_alpha_0(2.5)
    , alpha_beta_0(5)
    , beta_beta_0(2.5)
    , mu0(-20.0)
    , lambda0(0.01)
    , splicing_mu0(0.7)
    , splicing_alpha_alpha_0(1)
    , splicing_beta_alpha_0(1)
    , splicing_alpha_beta_0(1)
    , splicing_beta_beta_0(1)
{
    // organize transcripts by TSS and TTS
    transcript_ids.resize(ts.size());
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        pos_t tss_pos = t->strand == strand_pos ? t->min_start : t->max_end;
        Interval tss(t->seqname, tss_pos, tss_pos, t->strand);
        tss_group[tss].push_back(t->id);

        Interval tss_tts(t->seqname, t->min_start, t->max_end, t->strand);
        tss_tts_group[tss_tts].push_back(t->id);

        transcript_ids[t->id] = t->transcript_id;
        gene_ids.push_back(t->gene_id);
    }

    n_tss = tss_group.size();
    Logger::info("%u transcription start sites", n_tss);

    tss_index.resize(ts.size());
    tss_gene_ids.resize(tss_group.size());
    tss_transcript_ids.resize(tss_group.size());
    tss_tids.resize(tss_group.size());
    unsigned int tss_id = 0;
    BOOST_FOREACH (const TSMapItem& i, tss_group) {
        BOOST_FOREACH (const unsigned int& j, i.second) {
            tss_gene_ids[tss_id].insert(gene_ids[j]);
            tss_transcript_ids[tss_id].insert(transcript_ids[j]);
            tss_tids[tss_id].push_back(j);
            tss_index[j] = tss_id;
        }
        ++tss_id;
    }

    lambda_sampler = new LambdaSliceSampler(mu_dof, alpha, beta, mu,
                                            tss_usage, tss_base_precision,
                                            replicate_condition);
    alpha_sampler = new AlphaSliceSampler(alpha_alpha_0, beta_alpha_0, lambda);
    beta_sampler = new BetaSliceSampler(alpha_beta_0, beta_beta_0, lambda);
    rng = gsl_rng_alloc(gsl_rng_mt19937);
}


SwitchTest::~SwitchTest()
{
    delete lambda_sampler;
    delete alpha_sampler;
    gsl_rng_free(rng);
}


void SwitchTest::add_replicate(const char* condition_name, SampleDB& replicate_data)
{
    data[condition_name].push_back(&replicate_data);
}


void SwitchTest::run(unsigned int num_samples, const std::vector<double>& quantiles)
{
    m = 0;
    BOOST_FOREACH (CondMapItem& i, data) m += i.second.size();

    build_sample_matrix();
    sample_denorm.resize(m);

    const char* task_name = "Sampling";
    Logger::push_task(task_name, num_samples);

    mu.resize(data.size(), n_tss);
    std::fill(mu.data().begin(), mu.data().end(), mu0);

    alpha = alpha_alpha_0 / beta_alpha_0;
    beta  = alpha_beta_0 / beta_beta_0;

    lambda.resize(n_tss);
    std::fill(lambda.begin(), lambda.end(), alpha / beta);

    tss_usage.resize(m, n_tss);
    tss_base_precision.resize(m, n_tss);
    repl_sample_idx.resize(m);
    std::fill(repl_sample_idx.begin(), repl_sample_idx.end(), 0);

    tss_usage_norm_factor.resize(m);
    tss_usage_norm_factor_work.resize(n_tss);

    splicing.resize(m, ts.size());
    splicing_mu.resize(data.size(), ts.size());
    splicing_lambda.resize(n_tss);
    init_splicing_lambda_priors();

    std::vector<double> splicing_work;
    std::vector<double> splicing_mu;
    double splicing_prec;

    // TODO burnin
    for (unsigned int sample_num = 0; sample_num < num_samples; ++sample_num) {
        sample_replicate_transcript_abundance();
        sample_condition_tss_usage();
        sample_condition_splicing();
        mu_samples.push_back(mu);
        Logger::get_task(task_name).inc();
    }

    for (unsigned int cond1 = 0; cond1 < data.size() - 1; ++cond1) {
        for (unsigned cond2 = cond1 + 1; cond2 < data.size(); ++cond2) {
            size_t fnlen = condition_name[cond1].size() +
                           condition_name[cond2].size() +
                           strlen("transcription__.tsv") + 1;
            char* fn = new char[fnlen];
            snprintf(fn, fnlen, "transcription_%s_%s.tsv",
                     condition_name[cond1].c_str(),
                     condition_name[cond2].c_str());
            FILE* out = fopen(fn, "w");
            output_mu_samples(out, quantiles, cond1, cond2);
            fclose(out);
            delete [] fn;
        }
    }


    Logger::pop_task(task_name);
}


void SwitchTest::build_sample_matrix()
{
    const char* task_name = "Loading quantification data";
    Logger::push_task(task_name, m);

    if (data.size() == 0) {
        Logger::abort("No data to test.");
    }
    else if (data.size() == 1) {
        Logger::abort("Only one condition provided to test.");
    }

    condition_name.clear();
    condition_replicates.clear();
    replicate_condition.resize(m);

    std::map<TranscriptID, unsigned int> tids;
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        tids[t->transcript_id] = t->id;
    }

    effective_lengths.resize(m, tids.size());
    read_counts.resize(m);

    // for each replicate build a sample matrix in which rows are indexed by transcript
    // and columns by sample number.
    unsigned int cond_idx = 0;
    unsigned int repl_idx = 0;
    BOOST_FOREACH (CondMapItem& i, data) {
        if (i.second.size() == 0) {
            Logger::abort("Condition %s has no replicates.", i.first.c_str());
        }

        condition_name.push_back(i.first);
        condition_replicates.push_back(std::vector<unsigned int>());
        BOOST_FOREACH (SampleDB*& j, i.second) {
            read_counts[repl_idx] = j->get_num_reads();
            condition_replicates.back().push_back(repl_idx);
            replicate_condition[repl_idx] = cond_idx;
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
                effective_lengths(repl_idx, tid->second) = k->effective_length;
            }
            Logger::get_task(task_name).inc();
            ++repl_idx;
        }
        ++cond_idx;
    }

    Logger::pop_task(task_name);
}


void SwitchTest::init_splicing_lambda_priors()
{
    unsigned int max_group_size = 1;
    BOOST_FOREACH (TSMapItem& ts_item, tss_group) {
        max_group_size = std::max<unsigned int>(max_group_size, ts_item.second.size());
    }

    std::vector<unsigned int> group_size_counts(max_group_size + 1, 0);
    BOOST_FOREACH (TSMapItem& ts_item, tss_group) {
        group_size_counts[ts_item.second.size()]++;
    }

    max_group_size = 2;
    while (max_group_size < group_size_counts.size() &&
           group_size_counts[max_group_size] >=
               constants::min_tss_group_isoforms_conditioning) ++max_group_size;

    double prior_alpha_mean = splicing_alpha_alpha_0 / splicing_beta_alpha_0;
    double prior_beta_mean = splicing_alpha_beta_0 / splicing_beta_beta_0;
    splicing_alpha.resize(max_group_size, prior_alpha_mean);
    splicing_beta.resize(max_group_size, prior_beta_mean);

    double init_prec = exp(prior_alpha_mean / prior_beta_mean);
    std::fill(splicing_lambda.begin(), splicing_lambda.end(), init_prec);
}


void SwitchTest::compute_tss_usage_splicing(unsigned int i, unsigned int k)
{
    matrix_row<matrix<float> > tss_usage_row(tss_usage, i);
    matrix_row<matrix<float> > tss_base_precision_row(tss_base_precision, i);
    matrix_column<matrix<float> > samples_col(samples[i], k);

    std::fill(tss_usage_row.begin(), tss_usage_row.end(), 0.0);
    std::fill(tss_base_precision_row.begin(), tss_base_precision_row.end(), 0.0);
    for (unsigned int j = 0; j < ts.size(); ++j) {
        tss_usage_row[tss_index[j]] += constants::zero_eps + samples_col[j];

        double a = constants::zero_eps + samples_col[j];
        double b = effective_lengths(i, j);
        double c = read_counts[i];
        double l = a * b * c / sample_denorm[i];
        tss_base_precision_row[tss_index[j]] += l;
    }

    // upper quartile normalization between replicates
    std::copy(tss_usage_row.begin(), tss_usage_row.end(), tss_usage_norm_factor_work.begin());
    std::sort(tss_usage_norm_factor_work.begin(), tss_usage_norm_factor_work.end());
    tss_usage_norm_factor[i] =
        gsl_stats_quantile_from_sorted_data(&tss_usage_norm_factor_work.at(0), 1,
                                            n_tss, 0.75);

    double z = tss_usage_norm_factor[0] / tss_usage_norm_factor[i];

    matrix_row<matrix<float> > splicing_row(splicing, i);
    std::copy(samples_col.begin(), samples_col.end(), splicing_row.begin());
    for (unsigned int j = 0; j < ts.size(); ++j) {
        splicing_row[j] = (splicing_row[j] + constants::zero_eps) /
                tss_usage_row[tss_index[j]];
    }

    BOOST_FOREACH (float& x, tss_usage_row) {
        x = log2(z * x);
    }
}


#if 0
void SwitchTest::compute_splicing(unsigned int i, unsigned int k)
{
    matrix_row<matrix<float> > splicing_row(splicing, i);
    matrix_column<matrix<float> > samples_col(samples[i], k);

    assert(splicing_row.size() == samples_col.size());

    std::copy(samples_col.begin(), samples_col.end(), splicing_row.begin());

    // normalize to get abundance, given tss_usage
    for (unsigned int j = 0; j < ts.size(); ++j) {
        splicing_row[j] /= tss_usage(i, tss_index[j]);
    }
}
#endif


double SwitchTest::sample_log_likelihood(unsigned int i, unsigned int k)
{
    double p = 0.0;

    // tss usage log-probabilities
    compute_tss_usage_splicing(i, k);
    for (unsigned int j = 0; j < tss_usage.size2(); ++j)  {
        p += gaussian_lnpdf(mu(replicate_condition[i], j),
                            (tss_base_precision(i, j) * lambda[j]) /
                            (tss_base_precision(i, j) + lambda[j]),
                            tss_usage(i, j));
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
        // This is an approximation as it is not conditioning on mu or lambda.
        repl_sample_idx[i] = gsl_rng_uniform_int(rng, samples[i].size2());

        sample_denorm[i] = 0.0;
        for (unsigned int j = 0; j < ts.size(); ++j) {
            sample_denorm[i] +=
                effective_lengths(i, j) * samples[i](j, repl_sample_idx[i]);
        }

        compute_tss_usage_splicing(i, repl_sample_idx[i]);
    }
}


void SwitchTest::sample_condition_tss_usage()
{
    // sample mu
    for (unsigned int i = 0; i < data.size(); ++i) {
        for (unsigned int k = 0; k < tss_usage.size2(); ++k) {
            unsigned num_cond_replicates = condition_replicates[i].size();

            // sample from mu given lambda and transcript abundance
            double sample_mu = 0.0;
            BOOST_FOREACH (unsigned int& j, condition_replicates[i]) {
                sample_mu += tss_usage(j, k);
            }
            sample_mu /= (double) num_cond_replicates;


            double weighted_lambda_ik = 0.0;
            BOOST_FOREACH (unsigned int& j, condition_replicates[i]) {
                weighted_lambda_ik +=
                        (lambda[k] * tss_base_precision(j, k)) /
                        (lambda[k] + tss_base_precision(j, k));
            }

            double numer = sample_mu * weighted_lambda_ik + mu0 * lambda0;
            double denom = lambda0 + weighted_lambda_ik;
            double posterior_mu = numer / denom;
            double posterior_sigma = sqrt(1.0 / denom);

            // sampling from a student's t under the mixture interpretation
            double V = ran_scaled_inv_chisq(rng, mu_dof, sq(posterior_sigma));
            mu(i, k) = posterior_mu + gsl_ran_gaussian(rng, sqrt(V));
        }
    }

    // sample lambda
    for (unsigned int i = 0; i < n_tss; ++i) {
        lambda[i] = lambda_sampler->sample(i, lambda[i]);
    }

    // sample alpha
    alpha = alpha_sampler->sample(beta, alpha);

    double lambda_sum = std::accumulate(lambda.begin(), lambda.end(), 0.0);
    beta = gsl_ran_gamma(rng,
                         alpha_beta_0 + n_tss * alpha,
                         1.0 / (beta_beta_0 + lambda_sum));
}


void SwitchTest::sample_condition_splicing()
{
    // sample mu (given precision)
    // TODO

    // sample precision (given mu)
    // TODO
}


void SwitchTest::output_mu_samples(FILE* out, const std::vector<double>& quantiles,
                                   unsigned int cond1, unsigned int cond2)
{
    // header
    char name[100];
    fprintf(out, "gene_ids\ttranscript_ids");
    BOOST_FOREACH (const double& quantile, quantiles) {
        snprintf(name, 100, "%0.0f", 100.0 * quantile);
        fprintf(out, "\tcond1_low_%s\tcond1_high_%s\tcond2_low_%s\tcond2_high_%s\tlog2fc_low_%s\tlog2fc_high_%s",
                name, name, name, name, name, name);
    }
    fprintf(out, "\n");

    const unsigned int num_samples = mu_samples.size();
    std::vector<double> cond1_samples(num_samples), cond2_samples(num_samples),
                        log2fc_samples(num_samples);

    for (unsigned int i = 0; i < n_tss; ++i) {
        for (unsigned int j = 0; j < num_samples; ++j) {
            cond1_samples[j] = mu_samples[j](cond1, i);
            cond2_samples[j] = mu_samples[j](cond2, i);
            log2fc_samples[j] = cond1_samples[j] - cond2_samples[j];
        }

        std::sort(log2fc_samples.begin(), log2fc_samples.end());
        std::sort(cond1_samples.begin(), cond1_samples.end());
        std::sort(cond2_samples.begin(), cond2_samples.end());

        bool first = true;
        BOOST_FOREACH (const GeneID& gene_id, tss_gene_ids[i]) {
            if (!first) {
                fprintf(out, ",");
            }
            fprintf(out, "%s", gene_id.get().c_str());
            first = false;
        }
        fprintf(out, "\t");

        first = true;
        BOOST_FOREACH (const TranscriptID& transcript_id, tss_transcript_ids[i]) {
            if (!first) {
                fprintf(out, ",");
            }
            fprintf(out, "%s", transcript_id.get().c_str());
            first = false;
        }

        BOOST_FOREACH (const double& quantile, quantiles) {
            fprintf(out, "\t%e\t%e\t%e\t%e\t%e\t%e",
                    gsl_stats_quantile_from_sorted_data(&cond1_samples.at(0),
                                                        1, num_samples, quantile),
                    gsl_stats_quantile_from_sorted_data(&cond1_samples.at(0),
                                                        1, num_samples, 1.0 - quantile),
                    gsl_stats_quantile_from_sorted_data(&cond2_samples.at(0),
                                                        1, num_samples, quantile),
                    gsl_stats_quantile_from_sorted_data(&cond2_samples.at(0),
                                                        1, num_samples, 1.0 - quantile),
                    gsl_stats_quantile_from_sorted_data(&log2fc_samples.at(0),
                                                        1, num_samples, quantile),
                    gsl_stats_quantile_from_sorted_data(&log2fc_samples.at(0),
                                                        1, num_samples, 1.0 - quantile));
        }
        fprintf(out, "\n");
    }
}


