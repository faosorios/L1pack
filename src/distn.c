#include "distn.h"

/* dpqr-functions for Laplace distribution */
static double pdf_laplace(double, double, double, int);
static double cdf_laplace(double, double, double, int, int);
static double quantile_laplace(double, double, double, int, int);
static double rand_laplace(double, double);
/* ..end declarations */

static double pdf_laplace(double x, double location, double scale, int log_pdf)
{ /* density of the Laplace distribution */
  double y;

  y = fabs(x - location) / scale;

  return (log_pdf ?
    -0.5 * M_LN2 - log(scale) - M_SQRT2 * y :
    M_SQRT1_2 * exp(-M_SQRT2 * y) / scale);
}

void dlaplace(int *n, double *y, double *x, double *location, int *nloc, double *scale, int *nscale, int *give_log)
{ /* univariate Laplace random generation */
  int i, nobs = *n, na = *nloc, nb = *nscale, log_pdf = *give_log;

  for (i = 0; i < nobs; i++)
    y[i] = pdf_laplace(x[i], location[i % na], scale[i % nb], log_pdf);
}

static double cdf_laplace(double x, double location, double scale, int lower, int log_cdf)
{ /* distribution function of the Laplace distribution */
  double val;

  x = (x - location);
  val = .5 + .5 * sign(x) * (1. - exp(- M_SQRT2 * fabs(x) / scale));

  if (log_cdf) {
    return log(lower ? val : 1. - val);
  } else {
    return (lower ? val : 1. - val);
  }
}

void plaplace(int *n, double *y, double *x, double *location, int *nloc, double *scale, int *nscale, int *lower_tail, int *log_p)
{ /* univariate Laplace random generation */
  int i, nobs = *n, na = *nloc, nb = *nscale, lower = *lower_tail, log_cdf = *log_p;

  for (i = 0; i < nobs; i++)
    y[i] = cdf_laplace(x[i], location[i % na], scale[i % nb], lower, log_cdf);
}

static double quantile_laplace(double p, double location, double scale, int lower, int log_prob)
{ /* quantile function of the Laplace distribution */
  double q, val;

  if (scale == 0.0)
    return location;

  if (log_prob)
    p = exp(p);

  if (!lower)
    p = 1. - p;

  q = p - 0.5;

  if (q == 0.0)
    return location;

  if (q < 0.0)
    val = location + scale * M_SQRT1_2 * log1p(2. * q);
  else
    val = location - scale * M_SQRT1_2 * log(1. - 2. * q);

  return val;
}

void qlaplace(int *n, double *y, double *p, double *location, int *nloc, double *scale, int *nscale, int *lower_tail, int *log_p)
{ /* univariate Laplace random generation */
  int i, nobs = *n, na = *nloc, nb = *nscale, lower = *lower_tail, log_prob = *log_p;

  for (i = 0; i < nobs; i++)
    y[i] = quantile_laplace(p[i], location[i % na], scale[i % nb], lower, log_prob);
}

static double rand_laplace(double location, double scale)
{ /* random variates from the Laplace distribution */
  double u;

  u = unif_rand() - .5; /* U(-1/2,1/2) */

  if (scale == 0.)
    return location;
  else
    return location + M_SQRT1_2 * scale * sign(u) * log(1. - 2. * fabs(u));
}

void rlaplace(int *n, double *x, double *location, int *nloc, double *scale, int *nscale)
{ /* univariate Laplace random generation */
  int i, nobs = *n, na = *nloc, nb = *nscale;

  GetRNGstate();
  for (i = 0; i < nobs; i++)
    x[i] = rand_laplace(location[i % na], scale[i % nb]);
  PutRNGstate();
}
