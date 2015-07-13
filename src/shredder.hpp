
#ifndef ISOLATOR_SHREDDER_HPP
#define ISOLATOR_SHREDDER_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/uniform_01.hpp>

#include "common.hpp"
#include "nlopt/nlopt.h"


// A fast, generic, univariate slice sampler implementatation.
class Shredder
{
    public:
        Shredder(double lower_limit, double upper_limit, double tolerance);
        virtual ~Shredder();

        // change the tolerance
        void set_tolerance(double tolerance);

        // Generate a new sample given the current value x0.
        double sample(rng_t& rng, double x0);

        // Generate a maximum likelihood estimate
        double optimize(double x0);

        // These are kept as public fields to ease debugging and diagnostics.
        double x_min, x_max;

    protected:
        // Probability function. Return the log-probibility at x, as well as the
        // derivative d.
        virtual double f(double x, double& d) = 0;

        // Bounds on the parameter being sampled over
        double lower_limit, upper_limit, tolerance;

        boost::random::uniform_01<double> random_uniform_01;

    private:
        double find_slice_edge(double x0, double slice_height, double lp0,
                               double d0, int direction);

        nlopt_opt opt;

        friend double shredder_opt_objective(unsigned int _n, const double* _x,
                                             double* _grad, void* data);
};


// Some common distribution functions, with derivatives
//

#if 0
class PoissonLogPdf
{
    public:
        double f(float lambda, unsigned int k);
        double df_dlambda(float lambda, unsigned int k);
};
#endif


class NormalLogPdf
{
    public:
        double f(double mu, double sigma, const double* xs, size_t n);
        float f(float mu, float sigma, const float* xs, size_t n);

        double df_dx(double mu, double sigma, const double* xs, size_t n);

        double df_dmu(double mu, double sigma, const double* xs, size_t n);
        float df_dmu(float mu, float sigma, const float* xs, size_t n);

        double df_dsigma(double mu, double sigma, const double* xs, size_t n);
        float df_dsigma(float mu, float sigma, const float* xs, size_t n);

        double f(double mu, double sigma, double x);
        double df_dx(double mu, double sigma, double x);
};


class LogNormalLogPdf
{
    public:
        double f(double mu, double sigma, const double* xs, size_t n);
        double df_dx(double mu, double sigma, double x);
        double df_dmu(double mu, double sigma, const double* xs, size_t n);
        double df_dsigma(double mu, double sigma, const double* xs, size_t n);
};


#if 0
class StudentsTLogPdf
{
    public:
        double f(double nu, double mu, double sigma, const double* xs, size_t n);
        float f(float nu, float mu, float sigma, const float* xs, size_t n);
        double df_dx(double nu, double mu, double sigma, const double* xs, size_t n);
        double df_dmu(double nu, double mu, double sigma, const double* xs, size_t n);
        float df_dmu(float nu, float mu, float sigma, const float* xs, size_t n);
        double df_dsigma(double nu, double mu, double sigma, const double* xs, size_t n);
};
#endif


#if 0
class GammaLogPdf
{
    public:
        double f(double alpha, double beta, const double* xs, size_t n);
        float f(float alpha, float beta, const float* xs, size_t n);
        double df_dx(double alpha, double beta, const double* xs, size_t n);
        double df_dalpha(double alpha, double beta, const double* xs, size_t n);
        double df_dbeta(double alpha, double beta, const double* xs, size_t n);
        float df_dbeta(float alpha, float beta, const float* xs, size_t n);
};
#endif


#if 0
// Gamma parameterized by mean and shape
class AltGammaLogPdf
{
    public:
        double f(double mean, double shape, const double* xs, size_t n);
        float f(float mean, float shape, const float* xs, size_t n);
        double f(double mean, double shape, double x);
        float f(float mean, float shape, float x);
        double df_dx(double mean, double shape, const double* xs, size_t n);
        double df_dx(double mean, double shape, double x);
        double df_dmean(double mean, double shape, const double* xs, size_t n);
        float df_dmean(float mean, float shape, const float* xs, size_t n);
        double df_dshape(double shape, double mean, const double* xs, size_t n);
        float df_dshape(float shape, float mean, const float* xs, size_t n);
        float df_dshape(float mean, float shape, float x);
};
#endif


class InvGammaLogPdf
{
    public:
        double f(double alpha, double beta, const double* xs, size_t n);
        double df_dx(double alpha, double beta, const double* xs, size_t n);
        double df_dalpha(double alpha, double beta, const double* xs, size_t n);
        double df_dbeta(double alpha, double beta, const double* xs, size_t n);
};


class SqInvGammaLogPdf
{
    public:
        double f(double alpha, double beta, const double* xs, size_t n);
        double df_dx(double alpha, double beta, const double* xs, size_t n);
        double df_dalpha(double alpha, double beta, const double* xs, size_t n);
        double df_dbeta(double alpha, double beta, const double* xs, size_t n);
};


class BetaLogPdf
{
    public:
        double f(double alpha, double beta, double x);
        double df_dx(double alpha, double beta, double x);
        double df_dgamma(double gamma, double c, double x);
};


class DirichletLogPdf
{
    public:
        double f(double alpha,
                 const boost::numeric::ublas::matrix<double>* mean,
                 const boost::numeric::ublas::matrix<double>* data,
                 size_t n, size_t m);

        double df_dalpha(double alpha,
                         const boost::numeric::ublas::matrix<double>* mean,
                         const boost::numeric::ublas::matrix<double>* data,
                         size_t n, size_t m);
};


// 1-d logistic normal
class LogisticNormalLogPdf
{
    public:
        double f(double mu, double sigma, double x);
        double df_dx(double mu, double sigma, double x);

        // TODO: do these when we need them
        //double df_dmu(double x, double mu, double sigma);
        //double df_dsigma(double x, double mu ,double sigma);
};


#endif

