/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef ISOLATOR_DIRICHLET_H
#define ISOLATOR_DIRICHLET_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/**
 * \file
 * \brief Function for dealing with dirichlet distributions.
 */

/** Perform maximum-likelihood fitting of a dirichlet distribution.
 * \param n Dimension of the dirichlet distribution.
 * \param m Number of data points.
 * \param X A matrix of m points on a n-degree simplex.
 * \param mu The fit mean parameter.
 * \param prec The fit precision parameter.
 */
void fit_dirichlet(size_t n, size_t m, const double* X, double* mu, double* prec);


#ifdef __cplusplus
}
#endif

#endif

