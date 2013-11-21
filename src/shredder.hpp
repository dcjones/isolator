
#ifndef ISOLATOR_SHREDDER_HPP
#define ISOLATOR_SHREDDER_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/uniform_01.hpp>

#include "common.hpp"


// A fast, generic, univariate slice sampler implementatation.
class Shredder
{
    public:
        Shredder(double lower_limit, double upper_limit);
        virtual ~Shredder();

        // Generate a new sample given the current value x0.
        double sample(double x0);

    protected:
        // Probability function. Return the log-probibility at x, as well as the
        // derivative d.
        virtual double f(double x, double& d) = 0;

        // Bounds on the parameter being sampled over
        double lower_limit, upper_limit;

        rng_t rng;
        boost::random::uniform_01<double> random_uniform_01;

    private:
        double find_slice_edge(double x0, double slice_height, double lp0,
                               double d0, int direction);
};


// Some common distribution functions, with derivatives


class NormalLogPdf
{
    public:
        double f(double mu, double sigma, const double* xs, size_t n);
};


class StudentsTLogPdf
{
    public:
        double f(double nu, double mu, double sigma, const double* xs, size_t n);
        double df_dx(double nu, double mu, double sigma, const double* xs, size_t n);
        double df_dmu(double nu, double mu, double sigma, const double* xs, size_t n);
        double df_dsigma(double nu, double mu, double sigma, const double* xs, size_t n);
};


class GammaLogPdf
{
    public:
        double f(double alpha, double beta, const double* xs, size_t n);
        double df_dx(double alpha, double beta, const double* xs, size_t n);
        double df_dalpha(double alpha, double beta, const double* xs, size_t n);
        double df_dbeta(double alpha, double beta, const double* xs, size_t n);
};


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


#endif

