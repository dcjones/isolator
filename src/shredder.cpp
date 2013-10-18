
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>

#include "logger.hpp"
#include "shredder.hpp"


Shredder::Shredder(double lower_limit, double upper_limit)
    : lower_limit(lower_limit)
    , upper_limit(upper_limit)
{
}


Shredder::~Shredder()
{
}


double Shredder::sample(double x0)
{
    double d0;
    double lp0 = f(x0, d0);

    double slice_height = log(random_uniform_01(rng)) + lp0;

    double x_min = find_slice_edge(x0, slice_height, lp0, d0, -1);
    double x_max = find_slice_edge(x0, slice_height, lp0, d0,  1);

    double x;
    while (true) {
        x = x_min + (x_max - x_min) * random_uniform_01(rng);
        double d;
        double lp = f(x, d);

        if (lp >= slice_height) break;
        else if (x > x0) x_max = x0;
        else             x_min = x0;
    }

    return x;
}


double Shredder::find_slice_edge(double x0, double slice_height,
                                 double lp0, double d0, int direction)
{
    const double lp_eps = 1e-4;
    const double d_eps  = 1e-8;
    const double x_eps  = 1e-6;

    double lp = lp0 - slice_height;
    double d = d0;
    double x = x0;
    double x_bound = x0;

    while (fabs(lp) > lp_eps && fabs(lp/d) > d_eps) {
        double x1 = x - lp / d;

        // if we are very close to the boundry, and this iteration moves us past
        // the boundry, just give up.
        if (direction < 0 && fabs(x - lower_limit) <= x_eps && (x1 < x || lp > 0.0)) break;
        if (direction > 0 && fabs(x - upper_limit) <= x_eps && (x1 > x || lp > 0.0)) break;

        if (lp > 0) x_bound = x;

        // if we are moving in the wrong direction (i.e. toward the other root),
        // use bisection to correct course.
        if (direction < 0) {
            if (x1 > x_bound) {
                if (lp > 0) {
                    x = finite(lower_limit) ?
                            (lower_limit + x) / 2 : x - fabs(x - x1);
                }
                else {
                    x = (x_bound + x) / 2;
                }
            }
            else if (x1 < lower_limit) x = (lower_limit + x) / 2;
            else x = x1;
        }
        else {
            if (x1 < x_bound) {
                if (lp > 0) {
                    x = finite(upper_limit) ?
                            (upper_limit + x) / 2 : x + fabs(x - x1);
                }
                else {
                    x = (x_bound + x) / 2;
                }
            }
            else if (x1 > upper_limit) x = (upper_limit + x) / 2;
            else x = x1;
        }

        lp = f(x, d) - slice_height;
    }

    return x;
}


static double sq(double x)
{
    return x * x;
}


static double lbeta(double x, double y)
{
    return lgamma(x) + lgamma(y) - lgamma(x + y);
}


double NormalLogPdf::f(double mu, double sigma, const double* xs, size_t n)
{
    double part1 = n * (-log(2 * M_PI)/2 - log(sigma));
    double part2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part2 += sq(xs[i] - mu) / (2 * sq(sigma));
    }

    return part1 - part2;
}


double StudentsTLogPdf::f(double nu, double mu, double sigma, const double* xs, size_t n)
{
    double part1 =
        n * (lgamma((nu + 1) / 2) - lgamma(nu / 2) - log(sqrt(nu * M_PI) * sigma));

    double part2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part2 += log1p(sq((xs[i] - mu) / sigma) / nu);
    }

    return part1 - ((nu + 1) / 2) * part2;
}


double StudentsTLogPdf::df_dx(double nu, double mu, double sigma, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += (2 * (xs[i] - mu) / sq(sigma) / nu) / (1 + sq((xs[i] - mu) / sigma) / nu);
    }

    return -((nu + 1) / 2) * part;
}


double StudentsTLogPdf::df_dmu(double nu, double mu, double sigma, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += (2 * (xs[i] - mu) / sq(sigma) / nu) / (1 + sq((xs[i] - mu) / sigma) / nu);
    }

    return ((nu + 1) / 2) * part;
}


double StudentsTLogPdf::df_dsigma(double nu, double mu, double sigma, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += (2 * sq((xs[i] - mu) / sigma) / (nu * sigma)) /
                    (1 + sq((xs[i] - mu) / sigma) / nu);
    }

    return ((nu + 1) / 2) * part - n / sigma;
}



double InvGammaLogPdf::f(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += (alpha + 1) * log(xs[i]) + beta / xs[i];
    }

    return n * (alpha * log(beta) - lgamma(alpha)) - part;
}


double InvGammaLogPdf::df_dx(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += beta / sq(xs[i]) - (alpha + 1) / xs[i];
    }

    return part;
}


double InvGammaLogPdf::df_dalpha(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += log(xs[i]);
    }

    return n * (log(beta) - boost::math::digamma(alpha)) - part;
}


double InvGammaLogPdf::df_dbeta(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        part += 1 / xs[i];
    }

    return n * (alpha / beta) - part;
}


double SqInvGammaLogPdf::f(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double x = xs[i] * xs[i];
        part += (alpha + 1) * log(x) + beta / x;
    }

    return n * (alpha * log(beta) - lgamma(alpha)) - part;
}


double SqInvGammaLogPdf::df_dx(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double x = xs[i] * xs[i];
        part += beta / sq(x) - (alpha + 1) / x;
    }

    return part;
}


double SqInvGammaLogPdf::df_dalpha(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double x = xs[i] * xs[i];
        part += log(x);
    }

    return n * (log(beta) - boost::math::digamma(alpha)) - part;
}


double SqInvGammaLogPdf::df_dbeta(double alpha, double beta, const double* xs, size_t n)
{
    double part = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double x = xs[i] * xs[i];
        part += 1 / x;
    }

    return n * (alpha / beta) - part;
}


double BetaLogPdf::f(double alpha, double beta, double x)
{
    return (alpha - 1) * log(x) + (beta - 1) * log(1 - x) - lbeta(alpha, beta);
}


double BetaLogPdf::df_dx(double alpha, double beta, double x)
{
    return (alpha - 1) / x - (beta - 1) / (1 - x);
}
