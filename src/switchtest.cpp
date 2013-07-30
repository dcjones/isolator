
#include <set>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>
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
    return -LN_SQRT_TWO_PI + log(sqrt(lambda)) - sq(x - mu) * lambda / 2.0;
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
    SliceSampler(double lower_limit, double upper_limit);
    virtual ~SliceSampler();


protected:
    // Log-probability at value x.
    virtual double lpr(double x) = 0;

    // Generate a new sample given the current value x0.
    double sample(double x0);


private:
    double find_slice_edge(double x0, double slice_height, double init_step);
    double lower_limit, upper_limit;

    gsl_rng* rng;
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

    // XXX: optimization, not sampling
    double slice_height = log(gsl_rng_uniform(rng)) + y0;

    // find slice extents
    double step = 0.1 * x0;
    double x_min = find_slice_edge(x0, slice_height, -step);
    double x_max = find_slice_edge(x0, slice_height, +step);

    // sample
    double x = 0.0;
    double y;
    // while (true) {
        x = x_min + (x_max - x_min) * gsl_rng_uniform(rng);
        y = lpr(x);

        // if (y < slice_height) {
        //     if (x < x0) x_min = x;
        //     else        x_max = x;
        // }
        // else break;
    // }

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
    LambdaSliceSampler(double& alpha,
                       double& beta,
                       boost::numeric::ublas::matrix<float>& mu,
                       boost::numeric::ublas::matrix<float>& tss_usage,
                       std::vector<unsigned int>& replicate_condition)
        : SliceSampler(1e-4, INFINITY)
        , alpha(alpha)
        , beta(beta)
        , mu(mu)
        , tss_usage(tss_usage)
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
                    lambda,
                    tss_usage(i, tss_idx));
        }

        return p;
    }

private:
    unsigned tss_idx;

    double& alpha;
    double& beta;
    boost::numeric::ublas::matrix<float>& mu;
    boost::numeric::ublas::matrix<float>& tss_usage;
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


SwitchTest::SwitchTest(TranscriptSet& ts)
    : ts(ts)
    , alpha_alpha_0(5)
    , beta_alpha_0(2.5)
    , alpha_beta_0(5)
    , beta_beta_0(2.5)
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

        transcript_ids.push_back(t->transcript_id);
        gene_ids.push_back(t->gene_id);
    }

    n_tss = tss_group.size();
    Logger::info("%u transcription start sites", n_tss);

    tss_index.resize(ts.size());
    tss_gene_ids.resize(tss_group.size());
    tss_transcript_ids.resize(tss_group.size());
    unsigned int tss_id = 0;
    BOOST_FOREACH (const TSMapItem& i, tss_group) {
        BOOST_FOREACH (const unsigned int& j, i.second) {
            tss_gene_ids[tss_id].insert(gene_ids[j]);
            tss_transcript_ids[tss_id].insert(transcript_ids[j]);
            tss_index[j] = tss_id;
        }
        ++tss_id;
    }

    lambda_sampler = new LambdaSliceSampler(alpha, beta, mu,
                                            tss_usage, replicate_condition);
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

    const char* task_name = "Sampling";
    Logger::push_task(task_name, num_samples);

    mu.resize(data.size(), n_tss);
    std::fill(mu.data().begin(), mu.data().end(), mu0);

    // TODO: tinker with these. Move them to constants
    alpha_alpha_0 = 5.0;
    beta_alpha_0  = 2.5;
    alpha_beta_0  = 5.0;
    beta_beta_0   = 2.5;

    alpha = alpha_alpha_0 / beta_alpha_0;
    beta  = alpha_beta_0 / beta_beta_0;

    lambda.resize(n_tss);
    std::fill(lambda.begin(), lambda.end(), alpha / beta);

    tss_usage.resize(m, n_tss);
    repl_sample_idx.resize(m);
    std::fill(repl_sample_idx.begin(), repl_sample_idx.end(), 0);

    tss_usage_norm_factor.resize(m);
    tss_usage_norm_factor_work.resize(n_tss);

    FILE* samples_out = fopen("samples.tsv", "w");
    fprintf(samples_out, "sample_num\tcondition\ttss_id\tx\n");

    FILE* switchtest_out = fopen("switchtest.tsv", "w");
    fprintf(switchtest_out, "sample_num\tcondition\ttss_id\tmu\tlambda\n");

    // TODO burnin
    for (unsigned int sample_num = 0; sample_num < num_samples; ++sample_num) {
        sample_replicate_transcript_abundance();
        sample_condition_tss_usage();
        Logger::get_task(task_name).inc();
        Logger::info("alpha = %e, beta = %e", alpha, beta);

        // Debugging output
        unsigned int debug_n_tss = 100;
        for (unsigned int i = 0; i < data.size(); ++i) {
            for (unsigned int j = 0; j < debug_n_tss; ++j) {
                fprintf(switchtest_out, "%u\t%u\t%u\t%e\t%e\n",
                    sample_num, i, j, (double) mu(i, j), (double) lambda[j]);
            }
        }

        for (unsigned int i = 0; i < tss_usage.size1(); ++i) {
            for (unsigned int j = 0; j < debug_n_tss; ++j) {
                fprintf(samples_out, "%u\t%u\t%u\t%e\n",
                        sample_num, replicate_condition[i], j, (double) tss_usage(i, j));
            }
        }

        mu_samples.push_back(mu);
    }

    fclose(samples_out);
    fclose(switchtest_out);

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

    /* What is this going to return? A table?
     *
     * The relevent information is
     *  mu/lambda samples
     *
     * If we output another goddamn database, we need to figure out a way to
     * query it in a reasonable way.
     *
     */
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
            }
            Logger::get_task(task_name).inc();
            ++repl_idx;
        }
        ++cond_idx;
    }

    Logger::pop_task(task_name);
}


void SwitchTest::compute_tss_usage(unsigned int i, unsigned int k)
{
    matrix_row<matrix<float> > row(tss_usage, i);
    std::fill(row.begin(), row.end(), 0.0);
    for (unsigned int j = 0; j < ts.size(); ++j) {
        row[tss_index[j]] += samples[i](j, k);
    }

    BOOST_FOREACH (float& x, row) {
        x = std::max<float>(x, constants::zero_eps);
    }

    std::copy(row.begin(), row.end(), tss_usage_norm_factor_work.begin());
    std::sort(tss_usage_norm_factor_work.begin(), tss_usage_norm_factor_work.end());
    tss_usage_norm_factor[i] =
        gsl_stats_quantile_from_sorted_data(&tss_usage_norm_factor_work.at(0), 1,
                                            n_tss, 0.75);

    double z = tss_usage_norm_factor[0] / tss_usage_norm_factor[i];

    BOOST_FOREACH (float& x, row) {
        x = log2(z * x);
    }
}


double SwitchTest::sample_log_likelihood(unsigned int i, unsigned int k)
{
    double p = 0.0;

    // tss usage log-probabilities
    compute_tss_usage(i, k);
    for (unsigned int j = 0; j < tss_usage.size2(); ++j)  {
        p += gaussian_lnpdf(mu(replicate_condition[i], j),
                            lambda[j], tss_usage(i, j));
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
        compute_tss_usage(i, repl_sample_idx[i]);
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

            double weighted_lambda_ik = num_cond_replicates * lambda[k];
            double numer = sample_mu * weighted_lambda_ik + mu0 * lambda0;
            double denom = lambda0 + weighted_lambda_ik;
            double posterior_mu = numer / denom;
            double posterior_sigma = sqrt(1.0 / denom);
            mu(i, k) = posterior_mu + gsl_ran_gaussian(rng, posterior_sigma);
        }
    }

    // sample lambda
    for (unsigned int i = 0; i < n_tss; ++i) {
        lambda[i] = lambda_sampler->sample(i, lambda[i]);
    }


    // sample alpha
    alpha = alpha_sampler->sample(beta, alpha);

    double lambda_sum = 0.0;
    unsigned int num_nonzero_tss = 0;
    for (unsigned int j = 0; j < n_tss; ++j) {
        // If I allow larger values, beta shrinks to a tiny tiny number.
        // I would think it work the other way arround.
        if (lambda[j] < 1e3) {
            num_nonzero_tss++;
            lambda_sum += lambda[j];
        }

        // unsigned int i;
        // for (i = 0; i < m; ++i) {
        //     if (tss_usage(i, j) > -20.0) break;
        // }
        // if (i < m) {
        //     num_nonzero_tss++;
        //     lambda_sum += lambda[j];
        // }
    }

    Logger::info("num_nonzero_tss: %u", num_nonzero_tss);
    Logger::info("a_post = %e, b_post = %e",
                         (alpha_beta_0 + num_nonzero_tss * alpha),
                         (beta_beta_0 + lambda_sum));
    // Logger::info("E[beta] = %e",
    //                      (alpha_beta_0 + n_tss * alpha) /
    //                      (beta_beta_0 + lambda_sum));
    beta = gsl_ran_gamma(rng,
                         alpha_beta_0 + num_nonzero_tss * alpha,
                         1.0 / (beta_beta_0 + lambda_sum));
    // Fuck. This is correctly sampling but beta is still shrinking to 0. Why the fuck.

    // Ok, the problem now is some lambdas get very very large so that
    // likelihood goes no infinity for genes that are stuck at zero.

    // How do we deal with this? I think we need to figure

    // What I want to do is just compute these values for genes that are
    // expressed beyond some cutoff. But that isn't right because then we
    // could increase the likelihood just by forcing things to zero expression

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





