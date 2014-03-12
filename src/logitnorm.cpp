
#include <cmath>
#include <cstdio>
#include <boost/numeric/ublas/matrix.hpp>

#include "logitnorm.hpp"

static double sq(double x)
{
    return x * x;
}

// Compute log(exp(x) + exp(y)), avoiding overflow/underflow.
static double logaddexp(double x, double y)
{
    double u = x - y;
    if (u > 0.0) {
        return x + log1p(exp(-u));
    } else if (u <= 0.0) {
        return y + log1p(exp(u));
    } else {
        return x + y;
    }
}

double logit(double x)
{
    return log(x / (1 - x));
}


double invlogit(double x)
{
    double expx = exp(x);
    return expx / (1.0 + expx);
}


static double logitnorm_pdf(double mu, double sigma, double x)
{
    double ans = 1.0 / (sigma * sqrt(2*M_PI));
    ans *= exp(- sq(logit(x) - mu) / (2*sq(sigma)));
    ans *= 1 / (x * (1 - x));
    return ans;
}


static double logitnorm_logpdf(double mu, double sigma, double x)
{
    return -log(sigma) - log(sqrt(2*M_PI)) -
           sq(log(x / (1 - x)) - mu) / (2 * sq(sigma)) -
           log(x) - log(1-x);
}


double logitnorm_approx_mean(double mu, double sigma)
{
    const double step = 0.002;
    double accum = 0.0;
    double norm = 0.0;
    for (double x = step; x <= 1.0 - step; x += step) {
        double fx = logitnorm_pdf(mu, sigma, x);
        accum += x * fx;
        norm += fx;
    }

    return accum / norm;
}


static double logitnorm_is_bimodal(double mu, double sigma)
{
    const double step = 0.01;
    size_t hill_count = 0;
    bool increasing = false;
    double fx_prev = -INFINITY;
    for (double x = step; x < 1.0; x += step) {
        double fx = logitnorm_logpdf(mu, sigma, x);
        if (fx >= fx_prev) {
            if (!increasing) ++hill_count;
            increasing = true;
        }
        else increasing = false;
        fx_prev = fx;
    }

    //if (hill_count > 1) {
        //fprintf(stderr, "here\n");
    //}

    return hill_count > 1;


    // We need something a little more robust that this. If I just evaluated 

#if 0
    double a = logitnorm_pdf(mu, sigma, 0.01);
    double b = logitnorm_pdf(mu, sigma, 0.10);

    double d = logitnorm_pdf(mu, sigma, 0.90);
    double e = logitnorm_pdf(mu, sigma, 0.99);

    return a > b && e > d;
#endif
}


// Estimate logistic normal std. dev.
double logitnorm_approx_sd(double mu, double sigma)
{
    // disallow values of sigma that would result in a multimodal distribution.
    if (logitnorm_is_bimodal(mu, sigma)) {
        return INFINITY;
    }

    double exp_mu = exp(mu);
    double mean = exp_mu / (1 + exp_mu);

    const double step = 0.005;
    double accum = 0.0;
    double norm = 0.0;
    for (double x = step; x <= 1.0 - step; x += step) {
        //double fx = logitnorm_logpdf(mu, sigma, x);
        //accum = isnan(accum) ?
            //2 * log(fabs(x - mean)) + fx : logaddexp(accum, 2 * log(fabs(x - mean)) + fx);
        //norm = isnan(norm) ?
            //fx : logaddexp(norm, fx);

        double fx = logitnorm_pdf(mu, sigma, x);
        accum += sq(x - mean) * fx;
        norm += fx;
    }

    //return exp((accum - norm) / 2);
    return sqrt(accum / norm);
}


// Find the approxmiate sigma corresponding to a std. dev. of sd.
double estimate_logitnorm_sd_sigma(double mu, double sd)
{
    double a = 0.1;
    double b = 15.0;

    while (fabs(b - a) > 0.01) {
        double c = (a + b) / 2;
        double c_sd = logitnorm_approx_sd(mu, c);

        //if (fabs(sd - c_sd) < 0.001) {
            //break;
        //}

        if (sd > c_sd) a = c;
        else           b = c;
    }

    double c = (a + b) / 2;
    double c_sd = logitnorm_approx_sd(mu, c);
    if (fabs(mu) < 5.0 && c > 1.4 && sd <= 0.03) {
        fprintf(stderr, "here: %f (%f %f %f) %f %f %f\n", mu, a, c, b, fabs(sd - c_sd), c_sd,
                logitnorm_approx_sd(mu, 1.4));
    }

    return (a + b) / 2;
}


const double mu_min = -6;
const double mu_max = 6;
const double mu_step = 0.01;

const double sd_min  = 0.005;
const double sd_max  = 0.5;
const double sd_step = 0.005;

boost::numeric::ublas::matrix<float> sd_sigma_lookup;


void logitnorm_sd_sigma_init()
{
    size_t m = round((mu_max - mu_min) / mu_step);
    size_t n = round((sd_max - sd_min) / sd_step);
    sd_sigma_lookup.resize(m, n, false);

    for (size_t i = 0; i < m; ++i) {
        double mu = mu_min + i * mu_step;
        for (size_t j = 0; j < n; ++j) {
            double sd = sd_min + j * sd_step;
            sd_sigma_lookup(i, j) = estimate_logitnorm_sd_sigma(mu, sd);
        }
    }
}


double logitnorm_sd_sigma(double mu, double sd)
{
    // bilinear interpolation

    size_t m = sd_sigma_lookup.size1(),
           n = sd_sigma_lookup.size2();

    double i = std::max<double>(0.0,
               std::min<double>(m - 1, (mu - mu_min) / mu_step));
    double j = std::max<double>(0.0,
               std::min<double>(n - 1, (sd - sd_min) / sd_step));

    size_t i0 = floor(i),
           i1 = ceil(i),
           j0 = floor(j),
           j1 = ceil(j);

    if (i0 == i1) {
        if (i0 > 0)  --i0;
        else         ++i1;
    }

    if (j0 == j1) {
        if (j0 > 0)  --j0;
        else         ++j1;
    }

    double f00 = sd_sigma_lookup(i0, j0),
           f01 = sd_sigma_lookup(i0, j1),
           f10 = sd_sigma_lookup(i1, j0),
           f11 = sd_sigma_lookup(i1, j1);

    double i_span = (i1 - i0) * mu_step,
           j_span = (j1 - j0) * sd_step;

    double a = f00 * (i1 - i) * mu_step * (j1 - j) * sd_step +
               f10 * (i - i0) * mu_step * (j1 - j) * sd_step +
               f01 * (i1 - i) * mu_step * (j - j0) * sd_step +
               f11 * (i - i0) * mu_step * (j - j0) * sd_step;

    return a / (i_span * j_span);


#if 0

    ///////////////

    // TODO: this should use bilinear interpolation
    size_t i = std::min<size_t>(sd_sigma_lookup.size1() - 1,
                                (mu - mu_min) / mu_step);
    size_t j = std::min<size_t>(sd_sigma_lookup.size2() - 1,
                                (sd - sd_min) / sd_step);
    double ans = sd_sigma_lookup(i, j);
    return std::max<double>(sd_sigma_lookup(i, j), 0.001);
#endif
}


