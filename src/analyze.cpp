
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/thread.hpp>

#include "analyze.hpp"
#include "constants.hpp"
#include "shredder.hpp"

using namespace boost::numeric::ublas;
using namespace boost::accumulators;



static void assert_finite(double x)
{
    if (!boost::math::isfinite(x)) {
        Logger::abort("%f found where finite value expected.", x);
    }
}



class BetaDistributionSampler : public Shredder
{
    public:
        BetaDistributionSampler()
            : Shredder(1e-16, 1.0)
        {}

        // What is x here? I should really be sampling over a/b, right?
        double sample(double a0, double b0, double prec,
                      double a_prior, double b_prior,
                      const double* data, size_t n)
        {
            this->a0 = a0;
            this->b0 = b0;
            this->prec = prec;
            this->a_prior = a_prior;
            this->b_prior = b_prior;
            this->data = data;
            this->n = n;
            return Shredder::sample(a0 / (a0 + b0));
        }

    private:
        double a0, b0;
        double prec;
        double a_prior, b_prior;
        const double* data;
        size_t n;
        BetaLogPdf beta_logpdf;

    protected:
        double f(double x, double &d)
        {
            double fx = 0.0;
            d = 0.0;

            // prior
            fx += beta_logpdf.f(a_prior, b_prior, x);
            d += beta_logpdf.df_dx(a_prior, b_prior, x);

            // likelihood
            for (size_t i = 0; i < n; ++i) {
                fx += beta_logpdf.f(x * prec, (1.0 - x) * prec, data[i]);
                d += beta_logpdf.df_dgamma(x, prec, data[i]);
            }

            return fx;
        };
};


class NormalMuSampler
{
    public:
        NormalMuSampler()
        {
        }

        double sample(double sigma, const double* xs, size_t n,
                      double prior_mu, double prior_sigma)
        {
            double prior_var = prior_sigma * prior_sigma;
            double var = sigma * sigma;

            double part = 1/prior_var  + n/var;
            double posterior_mu =
                (prior_mu / prior_var + std::accumulate(xs, xs + n, 0.0) / var) / part;
            double posterior_sigma = sqrt(1 / part);

            return posterior_mu + random_normal(rng) * posterior_sigma;
        }

    private:
        rng_t rng;
        boost::random::normal_distribution<double> random_normal;
};


class NormalSigmaSampler
{
    public:
        NormalSigmaSampler()
        {
        }

        double sample(const double* xs, size_t n, double prior_alpha, double prior_beta)
        {
            double posterior_alpha = prior_alpha + n / 2.0;

            double part = 0.0;
            for (size_t i = 0; i < n; ++i) {
                part += xs[i] * xs[i];
            }
            double posterior_beta = prior_beta + part / 2.0;

            boost::random::gamma_distribution<double> dist(posterior_alpha, 1/posterior_beta);

            return sqrt(1 / dist(rng));
        }

    private:
        rng_t rng;
};


class TgroupMuSigmaSamplerThread
{
    public:
        TgroupMuSigmaSamplerThread(const matrix<double>& ts,
                                   matrix<double>& mu,
                                   const std::vector<double>& experiment_tgroup_mu,
                                   const std::vector<double>& experiment_tgroup_sigma,
                                   std::vector<double>& sigma,
                                   double& alpha,
                                   double& beta,
                                   const std::vector<int>& condition,
                                   const std::vector<std::vector<int> >& condition_samples,
                                   Queue<int>& tgroup_queue,
                                   Queue<int>& notify_queue)
            : ts(ts)
            , mu(mu)
            , experiment_tgroup_mu(experiment_tgroup_mu)
            , experiment_tgroup_sigma(experiment_tgroup_sigma)
            , sigma(sigma)
            , alpha(alpha)
            , beta(beta)
            , condition(condition)
            , condition_samples(condition_samples)
            , tgroup_queue(tgroup_queue)
            , notify_queue(notify_queue)
            , thread(NULL)
        {
            K = ts.size1();
            T = sigma.size();
            C = condition_samples.size();
            xs.resize(K);
        }

        void run()
        {
            int tgroup;
            while (true) {
                if ((tgroup = tgroup_queue.pop()) == -1) break;

                // sample mu
                for (size_t i = 0; i < C; ++i) {
                    size_t l = 0;
                    BOOST_FOREACH (int j, condition_samples[i]) {
                        xs[l++] = ts(j, tgroup);
                    }

                    mu(i, tgroup) = mu_sampler.sample(sigma[tgroup], &xs.at(0), l,
                                                      experiment_tgroup_mu[tgroup],
                                                      experiment_tgroup_sigma[tgroup]);
                    assert_finite(mu(i, tgroup));
                }

                // sample sigma
                for (size_t i = 0; i < K; ++i) {
                    xs[i] = ts(i, tgroup) - mu(condition[i], tgroup);
                }

                sigma[tgroup] = sigma_sampler.sample(&xs.at(0), K, alpha, beta);
                assert_finite(sigma[tgroup]);

                notify_queue.push(1);
            }
        }

        void start()
        {
            thread = new boost::thread(boost::bind(&TgroupMuSigmaSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        const matrix<double>& ts;
        matrix<double>& mu;
        const std::vector<double>& experiment_tgroup_mu;
        const std::vector<double>& experiment_tgroup_sigma;
        std::vector<double>& sigma;
        double& alpha;
        double& beta;
        const std::vector<int>& condition;
        const std::vector<std::vector<int> >& condition_samples;
        Queue<int>& tgroup_queue;
        Queue<int>& notify_queue;
        boost::thread* thread;

        // temporary data vector
        std::vector<double> xs;

        // number of replicates
        size_t K;

        // number of tgroups
        size_t T;

        // number of conditions
        size_t C;

        NormalMuSampler    mu_sampler;
        NormalSigmaSampler sigma_sampler;
};


// Sample parameters giving the mean within-group splicing proportions per
// condition.
class SpliceMuSigmaSamplerThread
{
    public:
        SpliceMuSigmaSamplerThread(
                        std::vector<std::vector<std::vector<double> > >& mu,
                        std::vector<std::vector<double> >& sigma,
                        const std::vector<std::vector<double> >& experiment_mu,
                        const std::vector<std::vector<double> >& experiment_sigma,
                        const double& splice_alpha,
                        const double& splice_beta,
                        const matrix<double>& Q,
                        const std::vector<unsigned int>& spliced_tgroup_indexes,
                        const std::vector<std::vector<unsigned int> >& tgroup_tids,
                        const std::vector<int>& condition,
                        const std::vector<std::vector<int> >& condition_samples,
                        Queue<int>& spliced_tgroup_queue,
                        Queue<int>& notify_queue)
            : mu(mu)
            , sigma(sigma)
            , experiment_mu(experiment_mu)
            , experiment_sigma(experiment_sigma)
            , splice_alpha(splice_alpha)
            , splice_beta(splice_beta)
            , Q(Q)
            , spliced_tgroup_indexes(spliced_tgroup_indexes)
            , tgroup_tids(tgroup_tids)
            , condition(condition)
            , condition_samples(condition_samples)
            , spliced_tgroup_queue(spliced_tgroup_queue)
            , notify_queue(notify_queue)
            , thread(NULL)
        {
            C = mu.size();
            K = Q.size1();
        }

        void run()
        {
            // temporary array for storing observation marginals
            std::vector<double> data(K);

            // temporary space for sampling precision
            size_t max_size2 = 0;
            BOOST_FOREACH (const std::vector<unsigned int>& tids, tgroup_tids) {
                max_size2 = std::max<size_t>(tids.size(), max_size2);
            }
            matrix<double> dataj(K, max_size2);
            matrix<double> _data_raw(K, max_size2);

            int j;
            while (true) {
                if ((j = spliced_tgroup_queue.pop()) == -1) break;
                size_t tgroup = spliced_tgroup_indexes[j];

                // transform quantification data for the selected tgroup
                // so we can treat it as normal (i.e. invert the logistic normal
                // transform)
                for (size_t i = 0; i < K; ++i) {
                    double datasum = 0.0;
                    for (size_t k = 0; k < tgroup_tids[tgroup].size(); ++k) {
                        unsigned int tid = tgroup_tids[tgroup][k];
                        dataj(i, k) = Q(i, tid);
                        _data_raw(i, k) = dataj(i, k);
                        datasum += dataj(i, k);
                    }

                    for (size_t k = 0; k < tgroup_tids[tgroup].size(); ++k) {
                        dataj(i, k) /= datasum;
                    }

                    size_t klast = tgroup_tids[tgroup].size() - 1;
                    for (size_t k = 0; k < tgroup_tids[tgroup].size(); ++k) {
                        dataj(i, k) = log(dataj(i, k) / dataj(i, klast));
                    }
                }

                // sample parameters for each condition
                for (size_t i = 0; i < C; ++i) {
                    for (size_t k = 0; k < tgroup_tids[tgroup].size() - 1; ++k) {
                        for (size_t l = 0; l < condition_samples[i].size(); ++l) {
                            size_t sample_idx = condition_samples[i][l];
                            data[l] = dataj(sample_idx, k);
                        }

                        mu[i][j][k] =
                            mu_sampler.sample(sigma[j][k],
                                              &data.at(0),
                                              condition_samples[i].size(),
                                              experiment_mu[j][k],
                                              experiment_sigma[j][k]);
                    }
                }

                for (size_t k = 0; k < tgroup_tids[tgroup].size() - 1; ++k) {
                    matrix_column<matrix<double> > col(dataj, k);
                    std::copy(col.begin(), col.end(), data.begin());
                    for (size_t l = 0; l < K; ++l) {
                        data[l] -= mu[condition[l]][j][k];
                    }

                    sigma[j][k] =
                        sigma_sampler.sample(&data.at(0), K,
                                             splice_alpha, splice_beta);
                }

                notify_queue.push(1);
            }
        }

        void start()
        {
            thread = new boost::thread(
                    boost::bind(&SpliceMuSigmaSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        // what we are sampling over
        std::vector<std::vector<std::vector<double> > >& mu;
        std::vector<std::vector<double> >& sigma;

        const std::vector<std::vector<double> >& experiment_mu;
        const std::vector<std::vector<double> >& experiment_sigma;

        const double& splice_alpha;
        const double& splice_beta;

        // current transcript quantification
        const matrix<double>& Q;

        const std::vector<unsigned int>& spliced_tgroup_indexes;
        const std::vector<std::vector<unsigned int> >& tgroup_tids;
        const std::vector<int>& condition;
        const std::vector<std::vector<int> >& condition_samples;

        Queue<int>& spliced_tgroup_queue;
        Queue<int>& notify_queue;

        size_t C, K;
        NormalMuSampler mu_sampler;
        NormalSigmaSampler sigma_sampler;
        rng_t rng;

        boost::thread* thread;
};


class ExperimentSpliceMuSigmaSamplerThread
{
    public:
        ExperimentSpliceMuSigmaSamplerThread(
                std::vector<std::vector<double> >& experiment_mu,
                std::vector<std::vector<double> >& experiment_sigma,
                const std::vector<std::vector<std::vector<double> > >&
                    condition_mu,
                const std::vector<unsigned int>& spliced_tgroup_indexes,
                const std::vector<std::vector<unsigned int> >& tgroup_tids,
                double experiment_prior_mu,
                double experiment_prior_sigma,
                double experiment_splice_alpha,
                double experiment_splice_beta,
                Queue<int>& spliced_tgroup_queue,
                Queue<int>& notify_queue)
            : experiment_mu(experiment_mu)
            , experiment_sigma(experiment_sigma)
            , condition_mu(condition_mu)
            , spliced_tgroup_indexes(spliced_tgroup_indexes)
            , tgroup_tids(tgroup_tids)
            , experiment_prior_mu(experiment_prior_mu)
            , experiment_prior_sigma(experiment_prior_sigma)
            , experiment_splice_alpha(experiment_splice_alpha)
            , experiment_splice_beta(experiment_splice_beta)
            , spliced_tgroup_queue(spliced_tgroup_queue)
            , notify_queue(notify_queue)
            , thread(NULL)
        {
            C = condition_mu.size();
        }

        void run()
        {
            // temporary space for marginals
            std::vector<double> data(C);

            // temporary space for sampling precision
            size_t max_size2 = 0;
            BOOST_FOREACH (const std::vector<unsigned int>& tids, tgroup_tids) {
                max_size2 = std::max<size_t>(tids.size(), max_size2);
            }
            matrix<double> meanj(C, max_size2);
            matrix<double> dataj(C, max_size2);

            int j;
            while (true) {
                if ((j = spliced_tgroup_queue.pop()) == -1) break;

                size_t tgroup = spliced_tgroup_indexes[j];

                for (size_t k = 0; k < tgroup_tids[tgroup].size() - 1; ++k) {
                    for (size_t i = 0; i < C; ++i) {
                        data[i] = condition_mu[i][j][k];
                    }

                    experiment_mu[j][k] =
                        mu_sampler.sample(experiment_sigma[j][k],
                                          &data.at(0), C,
                                          experiment_prior_mu,
                                          experiment_prior_sigma);

                    for (size_t i = 0; i < C; ++i) {
                        data[i] -= experiment_mu[j][k];
                    }

                    experiment_sigma[j][k] =
                        sigma_sampler.sample(&data.at(0), C,
                                             experiment_splice_alpha,
                                             experiment_splice_beta);
                }

                notify_queue.push(1);
            }
        }

        void start()
        {
            thread = new boost::thread(
                boost::bind(&ExperimentSpliceMuSigmaSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        std::vector<std::vector<double> >& experiment_mu;
        std::vector<std::vector<double> >& experiment_sigma;
        const std::vector<std::vector<std::vector<double> > >&
            condition_mu;
        const std::vector<unsigned int>& spliced_tgroup_indexes;
        const std::vector<std::vector<unsigned int> >& tgroup_tids;
        double experiment_prior_mu;
        double experiment_prior_sigma;
        double experiment_splice_alpha;
        double experiment_splice_beta;
        Queue<int>& spliced_tgroup_queue;
        Queue<int>& notify_queue;

        size_t C;
        NormalMuSampler mu_sampler;
        NormalSigmaSampler sigma_sampler;
        rng_t rng;

        boost::thread* thread;
};


class AlphaSampler : public Shredder
{
    public:
        AlphaSampler()
            : Shredder(1e-16, 1e5)
        {
        }

        double sample(double alpha0, double beta,
                      double alpha_alpha, double beta_alpha,
                      const double* sigmas, size_t n)
        {
            this->beta = beta;
            this->alpha_alpha = alpha_alpha;
            this->beta_alpha = beta_alpha;
            this->sigmas = sigmas;
            this->n = n;
            return Shredder::sample(alpha0);
        }

    private:
        double beta;
        double alpha_alpha;
        double beta_alpha;
        const double* sigmas;
        size_t n;

        InvGammaLogPdf prior_logpdf;
        SqInvGammaLogPdf likelihood_logpdf;

    protected:
        double f(double alpha, double &d)
        {
            d = 0.0;
            double fx = 0.0;

            d += prior_logpdf.df_dx(alpha_alpha, beta_alpha, &alpha, 1);
            fx += prior_logpdf.f(alpha_alpha, beta_alpha, &alpha, 1);

            d += likelihood_logpdf.df_dalpha(alpha, beta, sigmas, n);
            fx += likelihood_logpdf.f(alpha, beta, sigmas, n);

            return fx;
        }
};


class BetaSampler : public Shredder
{
    public:
        BetaSampler()
            : Shredder(1e-16, 1e5)
        {}

        double sample(double beta0, double alpha,
                      double alpha_beta, double beta_beta,
                      const double* sigmas, size_t n)
        {
            this->alpha = alpha;
            this->alpha_beta = alpha_beta;
            this->beta_beta = beta_beta;
            this->sigmas = sigmas;
            this->n = n;
            double ans = Shredder::sample(beta0);
            return ans;
            //return Shredder::sample(beta0);
        }

    private:
        double alpha;
        double alpha_beta;
        double beta_beta;
        const double* sigmas;
        size_t n;

        InvGammaLogPdf prior_logpdf;
        SqInvGammaLogPdf likelihood_logpdf;

    protected:
        double f(double beta, double &d)
        {
            d = 0.0;
            double fx = 0.0;

            d += prior_logpdf.df_dx(alpha_beta, beta_beta, &beta, 1);
            fx += prior_logpdf.f(alpha_beta, beta_beta, &beta, 1);

            d += likelihood_logpdf.df_dbeta(alpha, beta, sigmas, n);
            fx += likelihood_logpdf.f(alpha, beta, sigmas, n);

            return fx;
        }
};


class ExperimentTgroupMuSigmaSamplerThread
{
    public:
        ExperimentTgroupMuSigmaSamplerThread(
                std::vector<double>& experiment_tgroup_mu,
                double prior_mu,
                double prior_sigma,
                std::vector<double>& experiment_tgroup_sigma,
                double& experiment_tgroup_alpha,
                double& experiment_tgroup_beta,
                matrix<double>& tgroup_mu,
                Queue<int>& tgroup_queue,
                Queue<int>& notify_queue)
            : experiment_tgroup_mu(experiment_tgroup_mu)
            , prior_mu(prior_mu)
            , prior_sigma(prior_sigma)
            , experiment_tgroup_sigma(experiment_tgroup_sigma)
            , experiment_tgroup_alpha(experiment_tgroup_alpha)
            , experiment_tgroup_beta(experiment_tgroup_beta)
            , tgroup_mu(tgroup_mu)
            , tgroup_queue(tgroup_queue)
            , notify_queue(notify_queue)
            , thread(NULL)
        {
        }

        void run()
        {
            int tgroup;
            size_t C = tgroup_mu.size1();

            while (true) {
                if ((tgroup = tgroup_queue.pop()) == -1) break;

                double sigma = experiment_tgroup_sigma[tgroup];

                // sample experiment_tgroup_mu[tgroup]

                double prior_var = prior_sigma * prior_sigma;
                double var = sigma * sigma;

                matrix_column<matrix<double> > col(tgroup_mu, tgroup);

                double part = 1/prior_var + C/var;
                double posterior_mu =
                    (prior_mu / prior_var + std::accumulate(col.begin(), col.end(), 0.0) / var) / part;
                double posterior_sigma = sqrt(1 / part);

                experiment_tgroup_mu[tgroup] = posterior_mu + random_normal(rng) * posterior_sigma;

                // sample experiment_tgroup_sigma[tgroup]

                double mu = experiment_tgroup_mu[tgroup];
                double posterior_alpha = experiment_tgroup_alpha + C/2;
                part = 0.0;
                for (size_t i = 0; i < C; ++i) {
                    part += (col[i] - mu) * (col[i] - mu);
                }
                double posterior_beta = experiment_tgroup_beta + part/2;

                boost::random::gamma_distribution<double> dist(posterior_alpha, 1/posterior_beta);

                experiment_tgroup_sigma[tgroup] = sqrt(1 / dist(rng));

                notify_queue.push(1);
            }
        }

        void start()
        {
            thread = new boost::thread(boost::bind(&ExperimentTgroupMuSigmaSamplerThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        std::vector<double>& experiment_tgroup_mu;
        double prior_mu, prior_sigma;
        std::vector<double>& experiment_tgroup_sigma;
        double& experiment_tgroup_alpha;
        double& experiment_tgroup_beta;
        matrix<double>& tgroup_mu;

        Queue<int>& tgroup_queue;
        Queue<int>& notify_queue;
        boost::thread* thread;
        rng_t rng;
        boost::random::normal_distribution<double> random_normal;
};



Analyze::Analyze(size_t burnin,
                 size_t num_samples,
                 TranscriptSet& transcripts,
                 const char* genome_filename,
                 bool run_gc_correction)
    : burnin(burnin)
    , num_samples(num_samples)
    , transcripts(transcripts)
    , genome_filename(genome_filename)
    , run_gc_correction(run_gc_correction)
    , K(0)
    , C(0)
    , N(0)
    , T(0)
{
    N = transcripts.size();
    T = transcripts.num_tgroups();

    // TODO: constants (also maybe command line options eventually)
    tgroup_alpha_alpha = 1.0;
    tgroup_beta_alpha  = 5.0;
    tgroup_alpha_beta  = 1.0;
    tgroup_beta_beta   = 5.0;

    splice_alpha_alpha = 50.0;
    splice_beta_alpha  = 15.0;

    splice_alpha_beta  = 50.0;
    splice_beta_beta   = 15.0;

    experiment_tgroup_mu0 = -10;
    experiment_tgroup_sigma0 = 5.0;

    experiment_splice_mu0 = 0;
    experiment_splice_sigma0 = 0.1;

    tgroup_expr.resize(T);
    tgroup_row_data.resize(T);
    scale_work.resize(N);

    alpha_sampler = new AlphaSampler();
    beta_sampler = new BetaSampler();

    tgroup_tids = transcripts.tgroup_tids();

    for (size_t i = 0; i < tgroup_tids.size(); ++i) {
        if (tgroup_tids[i].size() > 1) {
            spliced_tgroup_indexes.push_back(i);
        }
    }

    Logger::info("Number of transcripts: %u", N);
    Logger::info("Number of transcription groups: %u", T);
    Logger::info("Number of tgroup with multiple isoforms: %u",
                  spliced_tgroup_indexes.size());

}


Analyze::~Analyze()
{
    delete alpha_sampler;
    delete beta_sampler;
}


void Analyze::add_sample(const char* condition_name, const char* filename)
{
    int c;
    std::map<std::string, int>::iterator it = condition_index.find(condition_name);
    if (it == condition_index.end()) {
        c = (int) condition_index.size();
        condition_index[condition_name] = c;
    }
    else c = (int) it->second;

    filenames.push_back(filename);
    condition.push_back(c);
    if (c >= (int) condition_samples.size()) condition_samples.resize(c + 1);
    condition_samples[c].push_back(K);
    ++K;
}


// Thread to initialize samplers and fragment models
class SamplerInitThread
{
    public:
        SamplerInitThread(const std::vector<std::string>& filenames, const char* fa_fn,
                          TranscriptSet& transcripts,
                          std::vector<FragmentModel*>& fms,
                          bool run_gc_correction,
                          std::vector<Sampler*>& samplers,
                          Queue<int>& indexes)
            : filenames(filenames)
            , fa_fn(fa_fn)
            , transcripts(transcripts)
            , fms(fms)
            , run_gc_correction(run_gc_correction)
            , samplers(samplers)
            , indexes(indexes)
            , thread(NULL)
        {
        }

        void run()
        {
            int index;
            while (true) {
                if ((index = indexes.pop()) == -1) break;

                fms[index] = new FragmentModel();
                fms[index]->estimate(transcripts, filenames[index].c_str(), fa_fn);

                samplers[index] = new Sampler(filenames[index].c_str(), fa_fn,
                                              transcripts, *fms[index], run_gc_correction);
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&SamplerInitThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        const std::vector<std::string>& filenames;
        const char* fa_fn;
        TranscriptSet& transcripts;
        std::vector<FragmentModel*>& fms;
        bool run_gc_correction;

        std::vector<Sampler*>& samplers;

        Queue<int>& indexes;

        boost::thread* thread;
};


// Threads to run sampler iterations
class SamplerTickThread
{
    public:
        SamplerTickThread(std::vector<Sampler*>& samplers,
                          matrix<double>& Q,
                          Queue<int>& tick_queue,
                          Queue<int>& tock_queue)
            : samplers(samplers)
            , Q(Q)
            , tick_queue(tick_queue)
            , tock_queue(tock_queue)
            , thread(NULL)
        { }

        void run()
        {
            int index;
            while (true) {
                if ((index = tick_queue.pop()) == -1) break;
                samplers[index]->sample();
                const std::vector<double>& state = samplers[index]->state();

                matrix_row<matrix<double> > row(Q, index);
                std::copy(state.begin(), state.end(), row.begin());

                // notify of completion
                tock_queue.push(1);
            }
        }

        void start()
        {
            if (thread != NULL) return;
            thread = new boost::thread(boost::bind(&SamplerTickThread::run, this));
        }

        void join()
        {
            thread->join();
            delete thread;
            thread = NULL;
        }

    private:
        std::vector<Sampler*> samplers;
        matrix<double>& Q;
        Queue<int>& tick_queue;
        Queue<int>& tock_queue;
        boost::thread* thread;

};


void Analyze::setup_samplers()
{
    fms.resize(K);
    qsamplers.resize(K);

    std::vector<SamplerInitThread*> threads(constants::num_threads);
    Queue<int> indexes;
    for (unsigned int i = 0; i < constants::num_threads; ++i) {
        threads[i] = new SamplerInitThread(filenames, genome_filename, transcripts,
                                           fms, run_gc_correction, qsamplers,
                                           indexes);
        threads[i]->start();
    }

    for (unsigned int i = 0; i < K; ++i) indexes.push(i);
    for (unsigned int i = 0; i < constants::num_threads; ++i) indexes.push(-1);

    for (unsigned int i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
        delete threads[i];
    }
}


void Analyze::setup_output(hid_t file_id)
{
    // transcript information
    // ----------------------
    {
        herr_t status;

        hsize_t dims[1] = { N };
        hid_t dataspace;
        if ((dataspace = H5Screate_simple(1, dims, NULL)) < 0) {
            Logger::abort("HDF5 dataspace creation failed.");
        }

        hid_t varstring_type = H5Tcopy(H5T_C_S1);
        if (varstring_type < 0 || H5Tset_size(varstring_type, H5T_VARIABLE) < 0) {
            Logger::abort("HDF5 type creation failed.");
        }

        const char** string_data = new const char* [N];

        // transcript_id table
        hid_t transcript_id_dataset =
            H5Dcreate2(file_id, "/transcript_id", varstring_type,
                       dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (transcript_id_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            string_data[t->id] = t->transcript_id.get().c_str();
        }

        status = H5Dwrite(transcript_id_dataset, varstring_type,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, string_data);
        if (status < 0) Logger::abort("HDF5 write failed.");

        H5Dclose(transcript_id_dataset);

        // gene_id table
        hid_t gene_id_dataset =
            H5Dcreate2(file_id, "/gene_id", varstring_type,
                       dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (gene_id_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            string_data[t->id] = t->gene_id.get().c_str();
        }

        status = H5Dwrite(gene_id_dataset, varstring_type,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, string_data);
        if (status < 0) Logger::abort("HDF5 write failed.");

        H5Dclose(gene_id_dataset);
        H5Tclose(varstring_type);

        delete [] string_data;

        // tgroup table
        hid_t tgroup_dataset =
            H5Dcreate2(file_id, "/tgroup", H5T_NATIVE_UINT,
                       dataspace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT);
        if (tgroup_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        unsigned int* tgroup_data = new unsigned int[N];

        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            tgroup_data[t->id] = t->tgroup;
        }

        status = H5Dwrite(tgroup_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                          H5S_ALL, tgroup_data);
        if (status < 0) {
            Logger::abort("HDF5 write failed.");
        }

        H5Dclose(tgroup_dataset);

        delete [] tgroup_data;
    }

    // sample quantification
    // ---------------------
    {
        herr_t status;

        hsize_t dims[3] = {num_samples, K, N};
        hsize_t chunk_dims[3] = {1, 1, N};

        hid_t dataset_create_property = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dataset_create_property, H5D_CHUNKED);
        H5Pset_chunk(dataset_create_property, 3, chunk_dims);
        H5Pset_deflate(dataset_create_property, 7);

        h5_sample_quant_dataspace = H5Screate_simple(3, dims, NULL);

        h5_sample_quant_dataset =
            H5Dcreate2(file_id, "/transcript_quantification", H5T_NATIVE_FLOAT,
                       h5_sample_quant_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);

        if (h5_sample_quant_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        H5Pclose(dataset_create_property);

        hsize_t sample_quant_mem_dims[2] = {K, N};
        h5_sample_quant_mem_dataspace =
            H5Screate_simple(2, sample_quant_mem_dims, NULL);

        hsize_t sample_quant_start[2] = {0, 0};
        status = H5Sselect_hyperslab(h5_sample_quant_dataspace, H5S_SELECT_SET,
                                     sample_quant_start, NULL,
                                     sample_quant_mem_dims, NULL);

        if (status < 0) {
            Logger::abort("HDF5 hyperslab creation failed.");
        }
    }

    // experiment parameters
    // ---------------------
    {
        herr_t status;

        if (H5Gcreate1(file_id, "/experiment", 0) < 0) {
            Logger::abort("HDF5 group creation failed.");
        }

        hsize_t dims[2] = {num_samples, T};
        hsize_t chunk_dims[2] = {1, T};

        hid_t dataset_create_property = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dataset_create_property, H5D_CHUNKED);
        H5Pset_chunk(dataset_create_property, 2, chunk_dims);
        H5Pset_deflate(dataset_create_property, 7);

        h5_experiment_tgroup_dataspace = H5Screate_simple(2, dims, NULL);

        h5_experiment_mean_dataset =
            H5Dcreate2(file_id, "/experiment/tgroup_mean", H5T_NATIVE_DOUBLE,
                       h5_experiment_tgroup_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);

        h5_experiment_sd_dataset =
            H5Dcreate2(file_id, "/experiment/tgroup_sd", H5T_NATIVE_DOUBLE,
                       h5_experiment_tgroup_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);

        if (h5_experiment_mean_dataset < 0 || h5_experiment_sd_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        hsize_t tgroup_row_dims = T;
        h5_tgroup_row_mem_dataspace = H5Screate_simple(1, &tgroup_row_dims, NULL);

        hsize_t tgroup_row_start = 0;
        status = H5Sselect_hyperslab(h5_tgroup_row_mem_dataspace, H5S_SELECT_SET,
                                     &tgroup_row_start, NULL, &tgroup_row_dims, NULL);

        // splicing parameters
        chunk_dims[1] = spliced_tgroup_indexes.size();
        H5Pset_chunk(dataset_create_property, 2, chunk_dims);

        h5_splice_param_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
        if (h5_splice_param_type < 0) {
            Logger::abort("HDF5 type creation failed.");
        }

        dims[1] = spliced_tgroup_indexes.size();
        h5_experiment_splice_dataspace = H5Screate_simple(2, dims, NULL);
        h5_splicing_mem_dataspace = H5Screate_simple(1, &dims[1], NULL);

        h5_experiment_splice_mu_dataset =
            H5Dcreate2(file_id, "/experiment/splice_mu", h5_splice_param_type,
                       h5_experiment_splice_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);
        if (h5_experiment_splice_mu_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        h5_experiment_splice_sigma_dataset =
            H5Dcreate2(file_id, "/experiment/splice_sigma", h5_splice_param_type,
                       h5_experiment_splice_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);
        if (h5_experiment_splice_sigma_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        H5Pclose(dataset_create_property);
    }

    // condition parameters
    // --------------------
    {
        if (H5Gcreate1(file_id, "/condition", 0) < 0) {
            Logger::abort("HDF5 group creation failed.");
        }

        hsize_t dims[3] = {num_samples, C, T};
        hsize_t chunk_dims[3] = {1, 1, T};

        hid_t dataset_create_property = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dataset_create_property, H5D_CHUNKED);
        H5Pset_chunk(dataset_create_property, 3, chunk_dims);
        H5Pset_deflate(dataset_create_property, 7);

        h5_condition_tgroup_dataspace = H5Screate_simple(3, dims, NULL);

        h5_condition_mean_dataset =
            H5Dcreate2(file_id, "/condition/tgroup_mean", H5T_NATIVE_FLOAT,
                       h5_condition_tgroup_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);

        if (h5_condition_mean_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        // splicing
        chunk_dims[2] = spliced_tgroup_indexes.size();
        H5Pset_chunk(dataset_create_property, 3, chunk_dims);

        dims[2] = spliced_tgroup_indexes.size();
        h5_condition_splice_mu_dataspace = H5Screate_simple(3, dims, NULL);

        h5_condition_splice_mu_dataset =
            H5Dcreate2(file_id, "/condition/splice_mu", h5_splice_param_type,
                       h5_condition_splice_mu_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);
        if (h5_condition_splice_mu_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        chunk_dims[1] = spliced_tgroup_indexes.size();
        H5Pset_chunk(dataset_create_property, 2, chunk_dims);

        dims[1] = spliced_tgroup_indexes.size();
        h5_condition_splice_sigma_dataspace = H5Screate_simple(2, dims, NULL);

        h5_condition_splice_sigma_dataset =
            H5Dcreate2(file_id, "/condition/splice_sigma", h5_splice_param_type,
                       h5_condition_splice_sigma_dataspace, H5P_DEFAULT,
                       dataset_create_property, H5P_DEFAULT);
        if (h5_condition_splice_sigma_dataset < 0) {
            Logger::abort("HDF5 dataset creation failed.");
        }

        H5Pclose(dataset_create_property);
    }

    h5_splice_work = new hvl_t[spliced_tgroup_indexes.size()];
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t num_tids = tgroup_tids[spliced_tgroup_indexes[i]].size();
        h5_splice_work[i].len = num_tids;
        h5_splice_work[i].p = reinterpret_cast<void*>(new float[num_tids]);
    }
}


void Analyze::cleanup()
{
    BOOST_FOREACH (FragmentModel* fm, fms) {
        delete fm;
    }
    fms.clear();

    BOOST_FOREACH (Sampler* qs, qsamplers) {
        delete qs;
    }
    qsamplers.clear();
}


void Analyze::qsampler_update_hyperparameters()
{
    for (size_t i = 0; i < K; ++i) {
        qsamplers[i]->hp.scale = scale[i];

        size_t c = condition[i];
        for (size_t j = 0; j < T; ++j) {
            qsamplers[i]->hp.tgroup_mu[j] = tgroup_mu(c, j);
            qsamplers[i]->hp.tgroup_sigma[j] = tgroup_sigma[j];
        }

        std::fill(qsamplers[i]->hp.splice_mu.begin(),
                  qsamplers[i]->hp.splice_mu.end(), 0.0);

        std::fill(qsamplers[i]->hp.splice_sigma.begin(),
                  qsamplers[i]->hp.splice_sigma.end(), 1.0);

        for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
            unsigned int tgroup = spliced_tgroup_indexes[j];
            for (size_t k = 0; k < tgroup_tids[tgroup].size(); ++k) {
                qsamplers[i]->hp.splice_mu[tgroup_tids[tgroup][k]] =
                    condition_splice_mu[c][j][k];
                qsamplers[i]->hp.splice_sigma[tgroup_tids[tgroup][k]] =
                    condition_splice_sigma[j][k];
            }
        }
    }
}


void Analyze::compute_ts()
{
    for (unsigned int i = 0; i < K; ++i) {
        matrix_row<matrix<double> > row(ts, i);
        std::fill(row.begin(), row.end(), 0.0);
        for (TranscriptSet::iterator t = transcripts.begin(); t != transcripts.end(); ++t) {
            row(t->tgroup) += Q(i, t->id);
        }

        for (size_t tgroup = 0; tgroup < T; ++tgroup) {
            row(tgroup) = log(row(tgroup));
            assert_finite(row(tgroup));
        }
    }
}


void Analyze::compute_xs()
{
    for (unsigned int i = 0; i < K; ++i) {
        std::fill(tgroup_expr.begin(), tgroup_expr.end(), 0.0);
        for (TranscriptSet::iterator t = transcripts.begin(); t != transcripts.end(); ++t) {
            tgroup_expr[t->tgroup] += Q(i, t->id);
        }

        for (TranscriptSet::iterator t = transcripts.begin(); t != transcripts.end(); ++t) {
            xs(i, t->id) = Q(i, t->id) / tgroup_expr[t->tgroup];
        }
    }
}


void Analyze::run(hid_t output_file_id)
{
    C = condition_index.size();
    Q.resize(K, N);
    ts.resize(K, T);
    xs.resize(K, N);
    scale.resize(K, 1.0);
    tgroup_sigma.resize(T);
    tgroup_mu.resize(K, T);
    experiment_tgroup_mu.resize(T);
    experiment_tgroup_sigma.resize(T);

    condition_splice_mu.resize(C);
    for (size_t i = 0; i < C; ++i) {
        condition_splice_mu[i].resize(spliced_tgroup_indexes.size());
        for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
            condition_splice_mu[i][j].resize(
                tgroup_tids[spliced_tgroup_indexes[j]].size());
            std::fill(condition_splice_mu[i][j].begin(),
                      condition_splice_mu[i][j].end(), 0.0);
        }
    }

    condition_splice_sigma.resize(spliced_tgroup_indexes.size());
    size_t flattened_sigma_size = 0;
    for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
        condition_splice_sigma[j].resize(
            tgroup_tids[spliced_tgroup_indexes[j]].size());
        std::fill(condition_splice_sigma[j].begin(),
                  condition_splice_sigma[j].end(), 1.0);
        flattened_sigma_size += condition_splice_sigma[j].size() - 1;
    }
    Logger::info("flattened_sigma_size = %lu", (unsigned long) flattened_sigma_size);

    condition_splice_sigma_work.resize(flattened_sigma_size);

    experiment_splice_mu.resize(spliced_tgroup_indexes.size());
    experiment_splice_sigma.resize(spliced_tgroup_indexes.size());
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        experiment_splice_mu[i].resize(
                tgroup_tids[spliced_tgroup_indexes[i]].size());
        std::fill(experiment_splice_mu[i].begin(),
                  experiment_splice_mu[i].end(), 0.0);
        experiment_splice_sigma[i].resize(
                tgroup_tids[spliced_tgroup_indexes[i]].size());
        std::fill(experiment_splice_sigma[i].begin(),
                  experiment_splice_sigma[i].end(), 0.1);
    }

    choose_initial_values();

    setup_samplers();
    setup_output(output_file_id);

    BOOST_FOREACH (Sampler* qsampler, qsamplers) {
        qsampler->start();
    }

    qsampler_threads.resize(constants::num_threads);
    BOOST_FOREACH (SamplerTickThread*& thread, qsampler_threads) {
        thread = new SamplerTickThread(qsamplers, Q, qsampler_tick_queue,
                                       qsampler_notify_queue);
        thread->start();
    }

    musigma_sampler_threads.resize(constants::num_threads);
    BOOST_FOREACH (TgroupMuSigmaSamplerThread*& thread, musigma_sampler_threads) {
        thread = new TgroupMuSigmaSamplerThread(ts, tgroup_mu,
                                                experiment_tgroup_mu, experiment_tgroup_sigma,
                                                tgroup_sigma,
                                                tgroup_alpha, tgroup_beta,
                                                condition, condition_samples,
                                                musigma_sampler_tick_queue,
                                                musigma_sampler_notify_queue);
        thread->start();
    }

    experiment_musigma_sampler_threads.resize(constants::num_threads);
    BOOST_FOREACH (ExperimentTgroupMuSigmaSamplerThread*& thread, experiment_musigma_sampler_threads) {
        thread = new ExperimentTgroupMuSigmaSamplerThread(
            experiment_tgroup_mu, experiment_tgroup_mu0, experiment_tgroup_sigma0,
            experiment_tgroup_sigma, experiment_tgroup_alpha, experiment_tgroup_beta,
            tgroup_mu, experiment_musigma_sampler_tick_queue,
            experiment_musigma_sampler_notify_queue);

        thread->start();
    }

    splice_mu_sigma_sampler_threads.resize(constants::num_threads);
    BOOST_FOREACH (SpliceMuSigmaSamplerThread*& thread, splice_mu_sigma_sampler_threads) {
        thread = new SpliceMuSigmaSamplerThread(
                condition_splice_mu, condition_splice_sigma,
                experiment_splice_mu, experiment_splice_sigma,
                splice_alpha, splice_beta, Q,
                spliced_tgroup_indexes,
                tgroup_tids,
                condition,
                condition_samples,
                splice_mu_sigma_sampler_tick_queue,
                splice_mu_sigma_sampler_notify_queue);
        thread->start();
    }

    experiment_splice_mu_sigma_sampler_threads.resize(constants::num_threads);
    BOOST_FOREACH (ExperimentSpliceMuSigmaSamplerThread*& thread,
                   experiment_splice_mu_sigma_sampler_threads) {
        thread = new ExperimentSpliceMuSigmaSamplerThread(
                experiment_splice_mu,
                experiment_splice_sigma,
                condition_splice_mu,
                spliced_tgroup_indexes,
                tgroup_tids,
                experiment_splice_mu0,
                experiment_splice_sigma0,
                experiment_splice_alpha,
                experiment_splice_beta,
                experiment_splice_mu_sigma_sampler_tick_queue,
                experiment_splice_mu_sigma_sampler_notify_queue);
        thread->start();
    }

    warmup();

    const char* task_name = "Sampling";
    Logger::push_task(task_name, burnin + num_samples);

    for (size_t i = 0; i < burnin; ++i) {
        sample();
        Logger::get_task(task_name).inc();
    }


    for (size_t i = 0; i < num_samples; ++i) {
        sample();
        write_output(i);
        Logger::get_task(task_name).inc();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        qsampler_tick_queue.push(-1);
        musigma_sampler_tick_queue.push(-1);
        experiment_musigma_sampler_tick_queue.push(-1);
        splice_mu_sigma_sampler_tick_queue.push(-1);
        experiment_splice_mu_sigma_sampler_tick_queue.push(-1);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        qsampler_threads[i]->join();
        musigma_sampler_threads[i]->join();
        experiment_musigma_sampler_threads[i]->join();
        splice_mu_sigma_sampler_threads[i]->join();
        experiment_splice_mu_sigma_sampler_threads[i]->join();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete qsampler_threads[i];
        delete musigma_sampler_threads[i];
        delete experiment_musigma_sampler_threads[i];
        delete splice_mu_sigma_sampler_threads[i];
        delete experiment_splice_mu_sigma_sampler_threads[i];
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        delete [] reinterpret_cast<float*>(h5_splice_work[i].p);
    }
    delete [] h5_splice_work;
    H5Dclose(h5_experiment_mean_dataset);
    H5Dclose(h5_experiment_sd_dataset);
    H5Sclose(h5_experiment_tgroup_dataspace);
    H5Dclose(h5_condition_mean_dataset);
    H5Sclose(h5_condition_tgroup_dataspace);
    H5Dclose(h5_sample_quant_dataset);
    H5Dclose(h5_experiment_splice_mu_dataset);
    H5Dclose(h5_experiment_splice_sigma_dataset);
    H5Dclose(h5_condition_splice_mu_dataset);
    H5Dclose(h5_condition_splice_sigma_dataset);
    H5Sclose(h5_sample_quant_dataspace);
    H5Sclose(h5_sample_quant_mem_dataspace);
    H5Sclose(h5_tgroup_row_mem_dataspace);
    H5Sclose(h5_experiment_splice_dataspace);
    H5Sclose(h5_condition_splice_mu_dataspace);
    H5Sclose(h5_condition_splice_sigma_dataspace);
    H5Sclose(h5_splicing_mem_dataspace);
    H5Tclose(h5_splice_param_type);

    Logger::pop_task(task_name);
    cleanup();
}


void Analyze::warmup()
{
    // An attempt at a more efficient warm-up procedure.

    // draw a few samples from the quantification samplers without priors
    for (size_t i = 0; i < 10; ++i) {
        for (size_t j = 0; j < K; ++j) {
            qsampler_tick_queue.push(j);
        }

        for (size_t j = 0; j < K; ++j) {
            qsampler_notify_queue.pop();
        }
    }

    BOOST_FOREACH (Sampler* sampler, qsamplers) {
        sampler->engage_priors();
    }

    compute_ts();
    compute_ts_scaling();
    compute_xs();

    // choose ml estimates for tgroup_mu
    for (size_t i = 0; i < C; ++i) {
        for (size_t j = 0; j < T; ++j) {
            double mu = 0.0;
            BOOST_FOREACH (int l, condition_samples[i]) {
                mu += ts(l, j);
            }
            tgroup_mu(i, j) = mu / condition_samples[i].size();
        }
    }

    // choose initially flat values for splicing parameters
    for (size_t i = 0; i < C; ++i) {
        for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
            std::fill(condition_splice_mu[i][j].begin(),
                      condition_splice_mu[i][j].end(), 0.0);
        }
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        std::fill(condition_splice_sigma[i].begin(),
                  condition_splice_sigma[i].end(), 1.0);
    }

    // initialially flat values for experiment splicing
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        std::fill(experiment_splice_mu[i].begin(),
                  experiment_splice_mu[i].end(), 0.0);

        std::fill(experiment_splice_sigma[i].begin(),
                  experiment_splice_sigma[i].end(), 0.1);
    }

    // ml estimates for experiment_tgroup_mu
    for (size_t j = 0; j < T; ++j) {
        double mu = 0.0;
        for (size_t i = 0; i < C; ++i) {
            mu += tgroup_mu(i, j);
        }
        experiment_tgroup_mu[j] = mu / C;
    }
}


void Analyze::sample()
{
    qsampler_update_hyperparameters();

    for (size_t i = 0; i < K; ++i) {
        qsampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < K; ++i) {
        qsampler_notify_queue.pop();
    }

    compute_ts_scaling();
    compute_ts();
    compute_xs();

    // sample condition-level parameters

    for (size_t i = 0; i < T; ++i) {
        musigma_sampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < T; ++i) {
        musigma_sampler_notify_queue.pop();
    }

    tgroup_alpha = alpha_sampler->sample(tgroup_alpha, tgroup_beta,
                                         tgroup_alpha_alpha, tgroup_beta_alpha,
                                         &tgroup_sigma.at(0), T);
    assert_finite(tgroup_alpha);

    tgroup_beta = beta_sampler->sample(tgroup_beta, tgroup_alpha,
                                       tgroup_alpha_beta, tgroup_beta_beta,
                                       &tgroup_sigma.at(0), T);
    assert_finite(tgroup_beta);

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        splice_mu_sigma_sampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        splice_mu_sigma_sampler_notify_queue.pop();
    }

    // sample experiment-level parameters

    for (size_t i = 0; i < T; ++i) {
        experiment_musigma_sampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < T; ++i) {
        experiment_musigma_sampler_notify_queue.pop();
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        experiment_splice_mu_sigma_sampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        experiment_splice_mu_sigma_sampler_notify_queue.pop();
    }

    experiment_tgroup_alpha =
        alpha_sampler->sample(experiment_tgroup_alpha,
                              experiment_tgroup_beta,
                              tgroup_alpha_alpha, tgroup_beta_alpha,
                              &experiment_tgroup_sigma.at(0), T);

    assert_finite(experiment_tgroup_alpha);

    experiment_tgroup_beta =
        beta_sampler->sample(experiment_tgroup_beta,
                             experiment_tgroup_alpha,
                             tgroup_alpha_beta, tgroup_beta_beta,
                             &experiment_tgroup_sigma.at(0), T);

    assert_finite(experiment_tgroup_beta);

    for (size_t i = 0, j = 0; j < condition_splice_sigma.size(); ++j) {
        for (size_t k = 0; k < condition_splice_sigma[j].size() - 1; ++k) {
            condition_splice_sigma_work[i++] = condition_splice_sigma[j][k];
        }
    }

    splice_alpha =
        alpha_sampler->sample(splice_alpha, splice_beta,
                              splice_alpha_alpha, splice_beta_alpha,
                              &condition_splice_sigma_work.at(0),
                              condition_splice_sigma_work.size());

    assert_finite(splice_alpha);

    /* I seem to get better results by fixing splice_beta and only sampling
     * splice_alpha. It's maybe a bit over-parameterized otherwise. */
#if 0
    splice_beta =
        beta_sampler->sample(splice_beta, splice_alpha,
                             splice_alpha_beta, splice_beta_beta,
                             &condition_splice_sigma_work.at(0),
                             condition_splice_sigma_work.size());
#endif
    splice_beta = 0.1;

    assert_finite(splice_beta);
}


void Analyze::write_output(size_t sample_num)
{
    hsize_t file_start2[2] = {sample_num, 0};
    hsize_t file_count2[2] = {1, T};
    herr_t status;

    status = H5Sselect_hyperslab(h5_experiment_tgroup_dataspace, H5S_SELECT_SET,
                                 file_start2, NULL, file_count2, NULL);

    if (status < 0) {
        Logger::abort("HD5 dataspace selection failed.");
    }

    std::copy(experiment_tgroup_mu.begin(), experiment_tgroup_mu.end(),
              tgroup_row_data.begin());
    status = H5Dwrite(h5_experiment_mean_dataset, H5T_NATIVE_FLOAT,
                      h5_tgroup_row_mem_dataspace, h5_experiment_tgroup_dataspace,
                      H5P_DEFAULT, &tgroup_row_data.at(0));

    if (status < 0) {
        Logger::abort("HD5 write operation failed.");
    }

    std::copy(experiment_tgroup_sigma.begin(), experiment_tgroup_sigma.end(),
              tgroup_row_data.begin());
    status = H5Dwrite(h5_experiment_sd_dataset, H5T_NATIVE_FLOAT,
                      h5_tgroup_row_mem_dataspace, h5_experiment_tgroup_dataspace,
                      H5P_DEFAULT, &tgroup_row_data.at(0));

    if (status < 0) {
        Logger::abort("HD5 write operation failed.");
    }


    hsize_t file_start3[3] = {sample_num, 0, 0};
    hsize_t file_count3[3] = {1, 1, T};
    status = H5Sselect_hyperslab(h5_condition_tgroup_dataspace, H5S_SELECT_SET,
                                 file_start2, NULL, file_count2, NULL);

    for (size_t i = 0; i < C; ++i) {
        file_start3[1] = i;
        status = H5Sselect_hyperslab(h5_condition_tgroup_dataspace, H5S_SELECT_SET,
                                     file_start3, NULL, file_count3, NULL);
        if (status < 0) {
            Logger::abort("HDF5 hyperslab selection failed.");
        }

        matrix_row<matrix<double > > mu_row(tgroup_mu, i);
        std::copy(mu_row.begin(), mu_row.end(), tgroup_row_data.begin());
        status = H5Dwrite(h5_condition_mean_dataset, H5T_NATIVE_FLOAT,
                          h5_tgroup_row_mem_dataspace, h5_condition_tgroup_dataspace,
                          H5P_DEFAULT, &tgroup_row_data.at(0));
        if (status < 0) {
            Logger::abort("HDF5 write operation failed.");
        }
    }


    hsize_t sample_quant_start[3] = {sample_num, 0, 0};
    hsize_t sample_quant_count[3] = {1, K, N};

    status = H5Sselect_hyperslab(h5_sample_quant_dataspace, H5S_SELECT_SET,
                                 sample_quant_start, NULL,
                                 sample_quant_count, NULL);

    if (status < 0) {
        Logger::abort("HDF5 dataspace selection failed.");
    }

    status = H5Dwrite(h5_sample_quant_dataset, H5T_NATIVE_DOUBLE,
                      h5_sample_quant_mem_dataspace, h5_sample_quant_dataspace,
                      H5P_DEFAULT, &Q.data()[0]);

    if (status < 0) {
        Logger::abort("HDF5 write operation failed.");
    }

    // write experiment and condition splicing parameters
    hsize_t experiment_splicing_start[2] = {sample_num, 0};
    hsize_t experiment_splicing_count[2] = {1, spliced_tgroup_indexes.size()};

    status = H5Sselect_hyperslab(h5_experiment_splice_dataspace,
                                 H5S_SELECT_SET,
                                 experiment_splicing_start, NULL,
                                 experiment_splicing_count, NULL);
    if (status < 0) {
        Logger::abort("HDF5 dataspace selection failed.");
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        float* xs = reinterpret_cast<float*>(h5_splice_work[i].p);
        for (size_t j = 0; j < h5_splice_work[i].len; ++j) {
            xs[j] = experiment_splice_mu[i][j];
        }
    }

    status = H5Dwrite(h5_experiment_splice_mu_dataset, h5_splice_param_type,
                      h5_splicing_mem_dataspace, h5_experiment_splice_dataspace,
                      H5P_DEFAULT, h5_splice_work);
    if (status < 0) {
        Logger::abort("HDF5 write operation failed.");
    }

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        float* xs = reinterpret_cast<float*>(h5_splice_work[i].p);
        for (size_t j = 0; j < h5_splice_work[i].len; ++j) {
            xs[j] = experiment_splice_sigma[i][j];
        }
    }

    status = H5Dwrite(h5_experiment_splice_sigma_dataset, h5_splice_param_type,
                      h5_splicing_mem_dataspace, h5_experiment_splice_dataspace,
                      H5P_DEFAULT, h5_splice_work);
    if (status < 0) {
        Logger::abort("HDF5 write operation failed.");
    }

    hsize_t condition_splice_mu_start[3] = {sample_num, 0, 0};
    hsize_t condition_splice_mu_count[3] = {1, 1, spliced_tgroup_indexes.size()};

    for (size_t i = 0; i < C; ++i) {
        for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
            float* xs = reinterpret_cast<float*>(h5_splice_work[j].p);
            for (size_t k = 0; k < h5_splice_work[j].len; ++k) {
                xs[k] = condition_splice_mu[i][j][k];
            }
        }

        condition_splice_mu_start[1] = i;
        status = H5Sselect_hyperslab(h5_condition_splice_mu_dataspace,
                                     H5S_SELECT_SET,
                                     condition_splice_mu_start, NULL,
                                     condition_splice_mu_count, NULL);
        if (status < 0) {
            Logger::abort("HDF5 dataspace selection failed.");
        }

        status = H5Dwrite(h5_condition_splice_mu_dataset, h5_splice_param_type,
                          h5_splicing_mem_dataspace, h5_condition_splice_mu_dataspace,
                          H5P_DEFAULT, h5_splice_work);

        if (status < 0) {
            Logger::abort("HDF5 write operation failed.");
        }
    }

    hsize_t condition_splice_sigma_start[2] = {sample_num, 0};
    hsize_t condition_splice_sigma_count[2] = {1, spliced_tgroup_indexes.size()};

    for (size_t j = 0; j < spliced_tgroup_indexes.size(); ++j) {
        float* xs = reinterpret_cast<float*>(h5_splice_work[j].p);
        for (size_t k = 0; k < h5_splice_work[j].len; ++k) {
            xs[k] = condition_splice_sigma[j][k];
        }
    }

    status = H5Sselect_hyperslab(h5_condition_splice_sigma_dataspace,
                                 H5S_SELECT_SET,
                                 condition_splice_sigma_start, NULL,
                                 condition_splice_sigma_count, NULL);
    if (status < 0) {
        Logger::abort("HDF5 dataspace selection failed.");
    }

    status = H5Dwrite(h5_condition_splice_sigma_dataset, h5_splice_param_type,
                      h5_splicing_mem_dataspace, h5_condition_splice_sigma_dataspace,
                      H5P_DEFAULT, h5_splice_work);
    if (status < 0) {
        Logger::abort("HDF5 write operation failed.");
    }
}


void Analyze::compute_ts_scaling()
{
    for (unsigned int i = 0; i < K; ++i) {
        matrix_row<matrix<double> > row(Q, i);

        // unscale abundance estimates so we can compute a new scale
        // and renormalize. I know this must seem weird.
        BOOST_FOREACH (double& x, row) x /= scale[i];

        // normalize according to an upper quantile
        std::copy(row.begin(), row.end(), scale_work.begin());
        std::sort(scale_work.begin(), scale_work.end());
        scale[i] = scale_work[scale_work.size() * 0.85];
    }

    for (int i = (int) K - 1; i >= 0; --i) {
        scale[i] = scale[0] / scale[i];
    }

    for (unsigned int i = 0; i < K; ++i) {
        matrix_row<matrix<double> > row(Q, i);
        BOOST_FOREACH (double& x, row) {
            x *= scale[i];
        }
    }
}


void Analyze::choose_initial_values()
{
    // tgroup_mu
    const double tgroup_mu_0 = -10;
    std::fill(tgroup_mu.data().begin(), tgroup_mu.data().end(), tgroup_mu_0);

    // tgroup_sigma
    const double tgroup_sigma_0 = 1.0;
    //const double tgroup_sigma_0 =
        //(tgroup_alpha_alpha / tgroup_beta_alpha) / (tgroup_alpha_beta / tgroup_beta_beta);
    std::fill(tgroup_sigma.begin(), tgroup_sigma.end(), tgroup_sigma_0);
    std::fill(experiment_tgroup_sigma.begin(), experiment_tgroup_sigma.end(), tgroup_sigma_0);

    tgroup_alpha = 0.1;
    tgroup_beta = 1.0;
    experiment_tgroup_alpha = 0.1;
    experiment_tgroup_beta = 1.0;

    splice_alpha = 2.0;
    splice_beta = 2.0;

    experiment_splice_alpha = 50;
    experiment_splice_beta = 15;
}

