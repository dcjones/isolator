
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/array.hpp>
#include <boost/thread.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>

#include "analyze.hpp"
#include "constants.hpp"
#include "shredder.hpp"

using namespace boost::numeric::ublas;
using namespace boost::accumulators;



static void assert_finite(double x)
{
    if (!finite(x)) {
        Logger::abort("%f found where finite value expected.", x);
    }
}


class TgroupMuSampler
{
    public:
        TgroupMuSampler()
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


class TgroupSigmaSampler
{
    public:
        TgroupSigmaSampler()
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
                                   double nu,
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
            , nu(nu)
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
        double nu;
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

        TgroupMuSampler    mu_sampler;
        TgroupSigmaSampler sigma_sampler;
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
            return Shredder::sample(beta0);
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
    tgroup_nu = 100.0;
    tgroup_alpha_alpha = 1.5;
    tgroup_beta_alpha  = 1.0;
    tgroup_alpha_beta  = 1.5;
    tgroup_beta_beta   = 1.0;

    experiment_tgroup_mu0 = -10;
    experiment_tgroup_sigma0 = 5.0;

    tgroup_expr.resize(T);
    scale_work.resize(N);

    alpha_sampler = new AlphaSampler();
    beta_sampler = new BetaSampler();

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
    herr_t status;
    hsize_t dims2[2] = {num_samples, 0};
    // hsize_t dims3[3] = {num_samples, 0, 0};

    if (H5Gcreate1(file_id, "/experiment", 0) < 0 ||
        H5Gcreate1(file_id, "/condition", 0) < 0 ||
        H5Gcreate1(file_id, "/sample", 0) < 0) {
        Logger::abort("Unable to create an HDF5 group.");
    }

    // data set creation properties
    const double fill_value = NAN;
    hid_t dataset_prop = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(dataset_prop, H5D_CONTIGUOUS);
    H5Pset_alloc_time(dataset_prop, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_value(dataset_prop, H5T_NATIVE_DOUBLE, &fill_value);
    H5Pset_fill_time(dataset_prop, H5D_FILL_TIME_ALLOC);

    // transcript data
    hsize_t transcript_dims = N;
    h5_transcript_dataspace = H5Screate_simple(1, &transcript_dims, NULL);

    hid_t h5_varstring_type = H5Tcopy(H5T_C_S1);
    if (h5_varstring_type < 0 || H5Tset_size(h5_varstring_type, H5T_VARIABLE) < 0) {
        Logger::abort("Failed to cerate HDF5 type.");
    }

    hid_t h5_transcript_id = H5Dcreate2(file_id, "/transcript_id", h5_varstring_type,
                                        h5_transcript_dataspace,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t h5_gene_id = H5Dcreate2(file_id, "/gene_id", h5_varstring_type,
                                  h5_transcript_dataspace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t h5_tgroup_id = H5Dcreate2(file_id, "/tgroup", H5T_NATIVE_UINT,
                                    h5_transcript_dataspace, H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);

    if (h5_transcript_id < 0 || h5_gene_id < 0 || h5_tgroup_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    // write transcript ids, gene ids
    {
        const char** transcript_string_data = new const char* [N];
        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            transcript_string_data[t->id] = t->transcript_id.get().c_str();
        }

        status = H5Dwrite(h5_transcript_id, h5_varstring_type,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, transcript_string_data);
        if (status < 0) Logger::abort("HDF5 write failed.");

        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            transcript_string_data[t->id] = t->gene_id.get().c_str();
        }

        status = H5Dwrite(h5_gene_id, h5_varstring_type,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, transcript_string_data);
        if (status < 0) Logger::abort("HDF5 write failed.");

        delete [] transcript_string_data;

        unsigned int* tgroup_data = new unsigned int[N];
        for (TranscriptSet::iterator t = transcripts.begin();
                t != transcripts.end(); ++t) {
            tgroup_data[t->id] = t->tgroup;
        }

        status = H5Dwrite(h5_tgroup_id, H5T_NATIVE_UINT,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, tgroup_data);
        if (status < 0) Logger::abort("HDF5 write failed.");

        delete[] tgroup_data;
    }

    H5Dclose(h5_transcript_id);
    H5Dclose(h5_gene_id);
    H5Dclose(h5_tgroup_id);
    H5Tclose(h5_varstring_type);

    // write tgroup ids


    // experiment parameters
    dims2[1] = T;
    h5_experiment_tgroup_dataspace = H5Screate_simple(2, dims2, NULL);

    h5_experiment_mean_dataset = H5Dcreate2(file_id, "/experiment/tgroup_mean",
                                            H5T_NATIVE_DOUBLE,
                                            h5_experiment_tgroup_dataspace,
                                            H5P_DEFAULT, dataset_prop, H5P_DEFAULT);

    h5_experiment_sd_dataset = H5Dcreate2(file_id, "/experiment/tgroup_sd",
                                          H5T_NATIVE_DOUBLE,
                                          h5_experiment_tgroup_dataspace,
                                          H5P_DEFAULT, dataset_prop, H5P_DEFAULT);

    if (h5_experiment_mean_dataset < 0 || h5_experiment_sd_dataset < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    hsize_t tgroup_row_dims = T;
    h5_tgroup_row_mem_dataspace = H5Screate_simple(1, &tgroup_row_dims, NULL);

    hsize_t tgroup_row_start = 0;
    status = H5Sselect_hyperslab(h5_tgroup_row_mem_dataspace, H5S_SELECT_SET,
                                 &tgroup_row_start, NULL, &tgroup_row_dims, NULL);

    // TODO: condition parameters

    // TODO: sample parameters
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
    experiment_tgroup_mu.resize(T);
    experiment_tgroup_sigma.resize(T);

    choose_initial_values();

    setup_samplers();

    const char* output_filename = "output.hdf5";
    hid_t output_file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (output_file_id < 0) {
        Logger::abort("Unable to open %s for writing.", output_filename);
    }
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
        thread = new TgroupMuSigmaSamplerThread(ts, tgroup_nu, tgroup_mu,
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
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        qsampler_threads[i]->join();
        musigma_sampler_threads[i]->join();
        experiment_musigma_sampler_threads[i]->join();
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete qsampler_threads[i];
        delete musigma_sampler_threads[i];
        delete experiment_musigma_sampler_threads[i];
    }

    H5Dclose(h5_experiment_mean_dataset);
    H5Dclose(h5_experiment_sd_dataset);
    H5Sclose(h5_experiment_tgroup_dataspace);
    H5Sclose(h5_tgroup_row_mem_dataspace);
    H5Fclose(output_file_id);

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

    // TODO: choose ml estimates for splicing paramaters


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

    // sample experiment-level parameters

    for (size_t i = 0; i < T; ++i) {
        experiment_musigma_sampler_tick_queue.push(i);
    }

    for (size_t i = 0; i < T; ++i) {
        experiment_musigma_sampler_notify_queue.pop();
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

    status = H5Dwrite(h5_experiment_mean_dataset, H5T_NATIVE_DOUBLE,
             h5_tgroup_row_mem_dataspace, h5_experiment_tgroup_dataspace,
             H5P_DEFAULT, &experiment_tgroup_mu.at(0));

    if (status < 0) {
        Logger::abort("HD5 write operation failed.");
    }

    H5Dwrite(h5_experiment_sd_dataset, H5T_NATIVE_DOUBLE,
             h5_tgroup_row_mem_dataspace, h5_experiment_tgroup_dataspace,
             H5P_DEFAULT, &experiment_tgroup_sigma.at(0));

    if (status < 0) {
        Logger::abort("HD5 write operation failed.");
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
}


