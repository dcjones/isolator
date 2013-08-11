
#include "dirichlet.h"
#include <nlopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

typedef struct
{
    size_t m;
    const double* X;
} fit_dirichlet_params;

static double dirichlet_ll_f(unsigned int n, const double* a,
                             double* grad, void* params)
{
    fit_dirichlet_params* p = (fit_dirichlet_params*) params;

    size_t i, j;

    /* likelihood */
    double L = 0.0;
    for (j = 0; j < p->m; ++j) {
        L += gsl_ran_dirichlet_lnpdf((size_t) n, a, p->X + j * n);
    }

    /* gradient */
    double a_sum = 0.0;
    for (i = 0; i < n; ++i) a_sum += a[i];
    double m_psi_a_sum = (double) p->m * (double) gsl_sf_psi(a_sum);

    for (i = 0; i < n; ++i) {
        grad[i] = m_psi_a_sum;
        grad[i] -= (double) p->m * gsl_sf_psi(a[i]);

        for (j = 0; j < p->m; ++j) {
            grad[i] += log(p->X[j * n + i]);
        }
    }

    return L;
}




void fit_dirichlet(size_t n, size_t m, const double* X, double* mu, double* prec)
{
    fit_dirichlet_params p;
    p.X = X;
    p.m = m;

    nlopt_opt fmax = nlopt_create(NLOPT_LD_LBFGS, n);
    nlopt_set_max_objective(fmax, dirichlet_ll_f, (void*) &p);
    nlopt_set_lower_bounds1(fmax, 1e-6);
    nlopt_set_upper_bounds1(fmax, 1e6);
    nlopt_set_initial_step1(fmax, 0.1);
    nlopt_set_ftol_rel(fmax, 1e-6);
    nlopt_set_maxeval(fmax, 2500);

    // initial value: simply sum all the datapoints
    size_t i, j;
    memset(mu, 0, n * sizeof(double));
    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
            mu[i] += X[j * n + i];
        }
    }

    double f_opt;
    nlopt_optimize(fmax, mu, &f_opt);

    /* reparameterize as mean and precision */
    *prec = 0;
    for (i = 0; i < n; ++i) *prec += mu[i];
    for (i = 0; i < n; ++i) mu[i] /= *prec;

    nlopt_destroy(fmax);
}


