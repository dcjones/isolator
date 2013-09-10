
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/thread.hpp>

#include "analyze.hpp"
#include "constants.hpp"
#include "shredder.hpp"

using namespace boost::numeric::ublas;



static void assert_finite(double x)
{
    if (!finite(x)) {
        Logger::abort("%f found where finite value expected.", x);
    }
}


class TgroupMuSampler : public Shredder
{
    public:
        TgroupMuSampler()
            : Shredder(-100, 100)
        {
        }

        double sample(double mu0, double nu, double sigma,
                      const std::vector<double>& xs, size_t n)
        {
            this->nu = nu;
            this->sigma = sigma;
            this->xs = &xs.at(0);
            this->n = n;
            return Shredder::sample(mu0);
        }

    private:
        double nu;
        double sigma;
        const double* xs;
        size_t n;

        StudentsTLogPdf logpdf;

    protected:
        double f(double mu, double &d)
        {
            // TODO: prior on mu
            d = logpdf.df_dmu(nu, mu, sigma, xs, n);
            return logpdf.f(nu, mu, sigma, xs, n);
        }
};


class TgroupSigmaSampler : public Shredder
{
    public:
        TgroupSigmaSampler()
            : Shredder(1e-16, 1e5)
        {
        }

        double sample(double sigma0, double nu,
                      const std::vector<double>& xs, double alpha, double beta,
                      size_t n)
        {
            this->nu = nu;
            this->xs = &xs.at(0);
            this->n = n;
            this->alpha = alpha;
            this->beta = beta;
            return Shredder::sample(sigma0);
        }

    private:
        double nu;
        const double* xs;
        size_t n;
        double alpha;
        double beta;

        StudentsTLogPdf likelihood_logpdf;
        InvGammaLogPdf prior_logpdf;

    protected:
        double f(double sigma, double &d)
        {
            d = 0.0;
            double fx = 0.0;

            d += prior_logpdf.df_dx(alpha, beta, &sigma, 1);
            fx += prior_logpdf.f(alpha, beta, &sigma, 1);

            d += likelihood_logpdf.df_dsigma(nu, 0.0, sigma, xs, n);
            fx += likelihood_logpdf.f(nu, 0.0, sigma, xs, n);

            return fx;
        }
};


class TgroupMuSigmaSamplerThread
{
    public:
        TgroupMuSigmaSamplerThread(const matrix<double>& ts,
                                   double nu,
                                   matrix<double>& mu,
                                   std::vector<double>& sigma,
                                   double& alpha,
                                   double& beta,
                                   const std::vector<int>& condition,
                                   const std::vector<std::vector<int> >& condition_samples,
                                   Queue<int>& tgroup_queue,
                                   Queue<int>& notify_queue)
            : ts(ts)
            , nu(nu)
            , mu(mu)
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

                    mu(i, tgroup) = mu_sampler.sample(mu(i, tgroup), nu,
                                                      sigma[tgroup], xs, l);
                    assert_finite(mu(i, tgroup));
                }

                // sample sigma
                for (size_t i = 0; i < K; ++i) {
                    xs[i] = ts(i, tgroup) - mu(condition[i], tgroup);
                }

                sigma[tgroup] = sigma_sampler.sample(sigma[tgroup], nu,
                                                     xs, alpha, beta, K);
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
        double nu;
        matrix<double>& mu;
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

        TgroupMuSampler    mu_sampler;
        TgroupSigmaSampler sigma_sampler;
};


class TgroupAlphaSampler : public Shredder
{
    public:
        TgroupAlphaSampler()
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
        InvGammaLogPdf likelihood_logpdf;

    protected:
        double f(double alpha, double &d)
        {
            d = 0.0;
            double fx = 0.0;

            d += prior_logpdf.df_dx(alpha_alpha, beta_alpha, &alpha, 1);
            fx += prior_logpdf.f(alpha_alpha, beta_alpha, &alpha, 1);

            d += likelihood_logpdf.df_dalpha(alpha, beta, sigmas, n);
            fx += prior_logpdf.f(alpha, beta, sigmas, n);

            return fx;
        }
};


class TgroupBetaSampler : public Shredder
{
    public:
        TgroupBetaSampler()
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
            return Shredder::sample(beta0);
        }

    private:
        double alpha;
        double alpha_beta;
        double beta_beta;
        const double* sigmas;
        size_t n;

        InvGammaLogPdf prior_logpdf;
        InvGammaLogPdf likelihood_logpdf;

    protected:
        double f(double beta, double &d)
        {
            d = 0.0;
            double fx = 0.0;

            d += prior_logpdf.df_dx(alpha_beta, beta_beta, &beta, 1);
            fx += prior_logpdf.f(alpha_beta, beta_beta, &beta, 1);

            d += likelihood_logpdf.df_dbeta(alpha, beta, sigmas, n);
            fx += prior_logpdf.f(alpha, beta, sigmas, n);

            return fx;
        }
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
    tgroup_nu = 100.0;
    tgroup_alpha_alpha = 5.0;
    tgroup_beta_alpha  = 2.5;
    tgroup_alpha_beta  = 5.0;
    tgroup_beta_beta   = 2.5;

    tgroup_expr.resize(T);

    alpha_sampler = new TgroupAlphaSampler();
    beta_sampler = new TgroupBetaSampler();

    Logger::info("Number of transcripts: %u", N);
    Logger::info("Number of transcription groups: %u", T);
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
                                              transcripts, *fms[index], run_gc_correction,
                                              true);
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
                samplers[index]->transition();
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


void Analyze::setup()
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
        qsamplers[i]->hp.tgroup_nu = tgroup_nu;

        size_t c = condition[i];
        for (size_t j = 0; j < T; ++j) {
            qsamplers[i]->hp.tgroup_mu[j] = tgroup_mu(c, j);
            qsamplers[i]->hp.tgroup_sigma[j] = tgroup_sigma[j];
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


void Analyze::run()
{
    C = condition_index.size();
    Q.resize(K, N);
    ts.resize(K, T);
    xs.resize(K, N);
    scale.resize(K, 1.0);
    tgroup_sigma.resize(T);
    tgroup_mu.resize(K, T);

    choose_initial_values();

    setup();
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
        thread = new TgroupMuSigmaSamplerThread(ts, tgroup_nu, tgroup_mu, tgroup_sigma,
                                                tgroup_alpha, tgroup_beta,
                                                condition, condition_samples,
                                                musigma_sampler_tick_queue,
                                                musigma_sampler_notify_queue);
        thread->start();
    }

    const char* task_name = "Sampling";
    Logger::push_task(task_name, burnin + num_samples);


    for (size_t i = 0; i < burnin; ++i) {
        sample();
        Logger::get_task(task_name).inc();
    }

    for (size_t i = 0; i < num_samples; ++i) {
        sample();
        // TODO: record the sample somehow
        Logger::get_task(task_name).inc();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        qsampler_tick_queue.push(-1);
        musigma_sampler_tick_queue.push(-1);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        qsampler_threads[i]->join();
        musigma_sampler_threads[i]->join();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete qsampler_threads[i];
        delete musigma_sampler_threads[i];
    }

    Logger::pop_task(task_name);
    cleanup();
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

    compute_ts();
    compute_xs();

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
}


#if 0
// Compute normalization constants for each sample.
void Analyze::compute_depth()
{
    const char* task_name = "Computing sample normalization constants";
    Logger::push_task(task_name, quantification.size() * (M/10));

    std::vector<double> sortedcol(N);
    depth.resize(K);

    // This is a posterior mean upper-quartile, normalized to the first sample
    // so depth tends to be close to 1.
    unsigned int i = 0;
    BOOST_FOREACH (matrix<float>& Q, quantification) {
        for (unsigned int j = 0; j < M; ++j) {
            matrix_column<matrix<float> > col(Q, i);
            std::copy(col.begin(), col.end(), sortedcol.begin());
            std::sort(sortedcol.begin(), sortedcol.end());
            depth[i] += gsl_stats_quantile_from_sorted_data(&sortedcol.at(0), 1, N, 0.75);
            if (j % 10 == 0) Logger::get_task(task_name).inc();
        }
        depth[i] /= M;
        ++i;
    }

    for (unsigned int i = 1; i < K; ++i) {
        depth[i] = depth[i] / depth[0];
    }
    depth[0] = 1.0;

    Logger::pop_task(task_name);
}
#endif


void Analyze::choose_initial_values()
{
    // tgroup_mu
    const double tgroup_mu_0 = -10;
    std::fill(tgroup_mu.data().begin(), tgroup_mu.data().end(), tgroup_mu_0);

    // tgroup_sigma
    const double tgroup_sigma_0 = 100.0;
    //const double tgroup_sigma_0 =
        //(tgroup_alpha_alpha / tgroup_beta_alpha) / (tgroup_alpha_beta / tgroup_beta_beta);
    std::fill(tgroup_sigma.begin(), tgroup_sigma.end(), tgroup_sigma_0);

    tgroup_alpha = 10.0;
    tgroup_beta = 2.0;
}

