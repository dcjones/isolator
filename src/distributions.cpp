
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "distributions.hpp"
#include "fastmath.hpp"


static float sq(float x)
{
    return x * x;
}


#if 0
static double cb(double x)
{
    return x * x * x;
}
#endif


void GammaLogPdf::set_alpha(float alpha)
{
    _alpha = alpha;
    lgamma_alpha = lgammaf(alpha);
}


void GammaLogPdf::set_beta(float beta)
{
    _beta = beta;
    log_beta = fastlog(beta);
}


float GammaLogPdf::f() const
{
    return _alpha * log_beta - lgamma_alpha - _beta * _x + (_alpha - 1.0f) * log_x;
}


float GammaLogPdf::x(float x)
{
    _x = x;
    log_x = fastlog(x);
    return f();
}


void GammaLogPdfDx::set_alpha(float alpha)
{
    _alpha = alpha;
}


void GammaLogPdfDx::set_beta(float beta)
{
    _beta = beta;
}


float GammaLogPdfDx::f() const
{
    return (_alpha - 1.0f) / _x - _beta;
}


float GammaLogPdfDx::x(float x)
{
    _x = x;
    return f();
}


void GammaLogPdfDBeta::set_alpha(float alpha)
{
    _alpha = alpha;
}


void GammaLogPdfDBeta::set_beta(float beta)
{
    _beta = beta;
}


float GammaLogPdfDBeta::f() const
{
    return _alpha / _beta - _x;
}


float GammaLogPdfDBeta::x(float x)
{
    _x = x;
    return f();
}


float AltGammaLogPdf::f() const
{
    return -(lgamma_shape + _shape * log_scale) + (_shape - 1.0f) * logx - _x / scale;
}


float AltGammaLogPdf::x(float x)
{
    _x = x;
    logx = fastlog(x);
    return f();
}


float AltGammaLogPdf::mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    return f();
}


void AltGammaLogPdf::set_mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
}


float AltGammaLogPdf::mean_x(float mean, float x)
{
    _mean = mean;
    _x = x;
    logx = fastlog(x);
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    return f();
}


float AltGammaLogPdf::shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    lgamma_shape = lgammaf(shape);
    log_scale = fastlog(scale);
    return f();
}


void AltGammaLogPdf::set_shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    lgamma_shape = lgammaf(shape);
    log_scale = fastlog(scale);
}


float AltGammaLogPdfDx::f() const
{
    return (_shape - 1.0f) / _x - 1.0f / scale;
}


float AltGammaLogPdfDx::x(float x)
{
    _x = x;
    return f();
}


float AltGammaLogPdfDx::mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
    return f();
}


void AltGammaLogPdfDx::set_mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
}


float AltGammaLogPdfDx::shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    return f();
}


void AltGammaLogPdfDx::set_shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
}


float AltGammaLogPdfDMean::f() const
{
    return _x * _shape / sq(_mean) - _shape / _mean;
}


float AltGammaLogPdfDMean::x(float x)
{
    _x = x;
    return f();
}


float AltGammaLogPdfDMean::mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
    return f();
}


float AltGammaLogPdfDMean::shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    return f();
}


float AltGammaLogPdfDShape::f() const
{
    return logx - _x / _mean - digamma_shape + log_scale * (_mean/sq(_shape));
}


float AltGammaLogPdfDShape::x(float x)
{
    _x = x;
    logx = fastlog(x);
    return f();
}


float AltGammaLogPdfDShape::mean(float mean)
{
    _mean = mean;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    return f();
}


float AltGammaLogPdfDShape::mean_x(float mean, float x)
{
    _x = x;
    logx = fastlog(x);
    _mean = mean;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    return f();
}


float AltGammaLogPdfDShape::shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    digamma_shape = boost::math::digamma(shape);
    return f();
}


void AltGammaLogPdfDShape::set_shape(float shape)
{
    _shape = shape;
    scale = _mean / _shape;
    log_scale = fastlog(scale);
    digamma_shape = boost::math::digamma(shape);
}


float PoissonLogPdf::f() const
{
    return _k * log_lambda - log_factorial_k - _lambda;
}


float PoissonLogPdf::k(unsigned int k)
{
    _k = k;
    log_factorial_k = lgammaf(k + 1);
    return f();
}


void PoissonLogPdf::set_k(unsigned int k)
{
    _k = k;
    log_factorial_k = lgammaf(k + 1);
}


float PoissonLogPdf::lambda(float lambda)
{
    _lambda = lambda;
    log_lambda = fastlog(lambda);
    return f();
}


float PoissonLogPdfDLambda::k(unsigned int k)
{
    _k = k;
    return _k / _lambda - 1.0f;
}


void PoissonLogPdfDLambda::set_k(unsigned int k)
{
    _k = k;
}


float PoissonLogPdfDLambda::lambda(float lambda)
{
    _lambda = lambda;
    return _k / _lambda - 1.0f;
}


void StudentsTLogPdf::set_nu(float nu)
{
    _nu = nu;
    nu_term = lgammaf((nu + 1.0f) / 2.0f) - lgammaf(nu / 2.0f);
    nu_sigma_term = - fastlog(sqrtf(nu * M_PI) * _sigma);
}


void StudentsTLogPdf::set_mu(float mu)
{
    _mu = mu;
}


void StudentsTLogPdf::set_sigma(float sigma)
{
    _sigma = sigma;
    nu_sigma_term = -fastlog(sqrtf(_nu * M_PI) * _sigma);
}


float StudentsTLogPdf::f() const
{
    return nu_term + nu_sigma_term -
        ((_nu + 1.0f) / 2.0f) * log1pf(sq((_x - _mu) / _sigma) / _nu);
}


float StudentsTLogPdf::x(float x)
{
    _x = x;
    return f();
}


void StudentsTLogPdfDx::set_nu(float nu)
{
    _nu = nu;
}


void StudentsTLogPdfDx::set_mu(float mu)
{
    _mu = mu;
}


void StudentsTLogPdfDx::set_sigma(float sigma)
{
    _sigma = sigma;
}


float StudentsTLogPdfDx::f() const
{
    float part = (2.0f * (_x - _mu) / sq(_sigma) / _nu) /
        (1.0f + sq((_x - _mu) / _sigma) / _nu);
    return -((_nu + 1.0f) / 2.0f) * part;
}


float StudentsTLogPdfDx::x(float x)
{
    _x = x;
    return f();
}


void StudentsTLogPdfDMu::set_nu(float nu)
{
    _nu = nu;
}


void StudentsTLogPdfDMu::set_mu(float mu)
{
    _mu = mu;
}


void StudentsTLogPdfDMu::set_sigma(float sigma)
{
    _sigma = sigma;
}


float StudentsTLogPdfDMu::f() const
{
    float part = (2.0f * (_x - _mu) / sq(_sigma) / _nu) /
        (1.0f + sq((_x - _mu) / _sigma) / _nu);
    return ((_nu + 1.0f) / 2.0f) * part;
}


float StudentsTLogPdfDMu::x(float x)
{
    _x = x;
    return f();
}
