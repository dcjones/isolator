
#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP


class GammaLogPdf
{
    public:
        GammaLogPdf()
            : _alpha(NAN), _beta(NAN), _x(NAN),
              lgamma_alpha(NAN), log_beta(NAN), log_x(NAN) {}

        void set_alpha(float alpha);
        void set_beta(float beta);
        float x(float x);

    private:
        float f() const;

        float _alpha, _beta, _x;
        float lgamma_alpha, log_beta, log_x;
};


class GammaLogPdfDx
{
    public:
        GammaLogPdfDx()
            : _alpha(NAN), _beta(NAN), _x(NAN) {}

        void set_alpha(float alpha);
        void set_beta(float beta);
        float x(float x);

    private:
        float f() const;

        float _alpha, _beta, _x;
};


class GammaLogPdfDBeta
{
    public:
        GammaLogPdfDBeta()
            : _alpha(NAN), _beta(NAN), _x(NAN) {}

        void set_alpha(float alpha);
        void set_beta(float beta);
        float x(float x);

    private:
        float f() const;

        float _alpha, _beta, _x;
};


// Gamma distribution parameterized by mean and shape
class AltGammaLogPdf
{
    public:
        AltGammaLogPdf()
            : _x(NAN), logx(NAN), _mean(NAN), _shape(NAN),
              lgamma_shape(NAN), scale(NAN), log_scale(NAN) {}

        float x(float x);
        float mean(float mean);
        float shape(float shape);
        float mean_x(float mean, float x);

        void set_mean(float mean);
        void set_shape(float shape);

    private:
        float f() const;

        float _x,
              logx,
              _mean,
              _shape,
              lgamma_shape,
              scale,
              log_scale;
};


class AltGammaLogPdfDx
{
    public:
        AltGammaLogPdfDx()
            : _x(NAN), _mean(NAN), _shape(NAN), scale(NAN) {}

        float x(float x);
        float mean(float mean);
        float shape(float shape);
        void set_mean(float mean);
        void set_shape(float shape);

    private:
        float f() const;
        float _x,
              _mean,
              _shape,
              scale;
};


class AltGammaLogPdfDMean
{
    public:
        AltGammaLogPdfDMean()
            : _x(NAN), _mean(NAN), _shape(NAN), scale(NAN) {}

        float x(float x);
        float mean(float mean);
        float shape(float shape);

    private:
        float f() const;
        float _x,
              _mean,
              _shape,
              scale;
};


class AltGammaLogPdfDShape
{
    public:
        AltGammaLogPdfDShape()
            : _x(NAN), logx(NAN), _mean(NAN), _shape(NAN),
              scale(NAN), log_scale(NAN), digamma_shape(NAN) {}

        float x(float x);
        float mean(float mean);
        float shape(float shape);
        void set_shape(float shape);
        float mean_x(float mean, float x);

    private:
        float f() const;

        float _x,
              logx,
              _mean,
              _shape,
              scale,
              log_scale,
              digamma_shape;
};


class PoissonLogPdf
{
    public:
        PoissonLogPdf()
            : _lambda(NAN), log_lambda(NAN), log_factorial_k(NAN) {}

        float k(unsigned int k);
        float lambda(float lambda);
        void set_k(unsigned int k);

    private:
        float f() const;

        unsigned int _k;
        float _lambda,
              log_lambda,
              log_factorial_k;
};


class PoissonLogPdfDLambda
{
    public:
        PoissonLogPdfDLambda()
            : _k(NAN), _lambda(NAN) {}

        float k(unsigned int k);
        float lambda(float lambda);
        void set_k(unsigned int k);

    private:
        float _k, _lambda;
};


class StudentsTLogPdf
{
    public:
        StudentsTLogPdf()
            : _x(NAN), _nu(NAN), _mu(NAN), _sigma(NAN),
              nu_term(NAN), nu_sigma_term(NAN) {}

        float x(float x);
        void set_nu(float nu);
        void set_mu(float mu);
        void set_sigma(float sigma);

    private:
        float f() const;

        float _x, _nu, _mu, _sigma;
        float nu_term,
              nu_sigma_term;
};


class StudentsTLogPdfDx
{
    public:
        StudentsTLogPdfDx()
            : _x(NAN), _nu(NAN), _mu(NAN), _sigma(NAN) {}

        float x(float x);
        void set_nu(float nu);
        void set_mu(float mu);
        void set_sigma(float sigma);

    private:
        float f() const;

        float _x, _nu, _mu, _sigma;
};


class StudentsTLogPdfDMu
{
    public:
        StudentsTLogPdfDMu()
            : _x(NAN), _nu(NAN), _mu(NAN), _sigma(NAN) {}

        float x(float x);
        void set_nu(float nu);
        void set_mu(float mu);
        void set_sigma(float sigma);

    private:
        float f() const;

        float _x, _nu, _mu, _sigma;
};


#endif


