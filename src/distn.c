/* ID: distn.c, last updated 2023-12-28, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

/* dpq-functions for Laplace distribution */
static double pdf_laplace(double, double, double, int);
static double cdf_laplace(double, double, double, int, int);
static double quantile_laplace(double, double, double, int, int);
/* ..end declarations */

/* ========================================================================== *
 * pdf, cdf and quantile functions for the univariate Laplace distribution
 * ========================================================================== */

static double pdf_laplace(double x, double location, double scale, int log_pdf)
{ /* density of the Laplace distribution */
  double y;

  y = fabs(x - location) / scale;

  return (log_pdf ?
    -.5 * M_LN2 - log(scale) - M_SQRT2 * y :
    M_SQRT1_2 * exp(-M_SQRT2 * y) / scale);
}

static double cdf_laplace(double x, double location, double scale, int lower, int log_cdf)
{ /* distribution function of the Laplace distribution */
  double val;

  x = (x - location);

  if (x < 0.)
    val = .5 * exp(M_SQRT2 * x / scale);
  else 
    val = 1. - .5 * exp(-M_SQRT2 * x / scale);

  if (log_cdf) {
    return log(lower ? val : 1. - val);
  } else {
    return (lower ? val : 1. - val);
  }
}

static double quantile_laplace(double p, double location, double scale, int lower, int log_prob)
{ /* quantile function of the Laplace distribution */
  double q, val;

  if (scale == 0.)
    return location;

  if (log_prob)
    p = exp(p);

  if (!lower)
    p = 1. - p;

  q = p;

  if (q == .5)
    return location;

  if (q < .5)
    val = location + scale * M_SQRT1_2 * log(2. * q);
  else
    val = location - scale * M_SQRT1_2 * log(2. * (1. - q));

  return val;
}

void d_laplace(int *n, double *y, double *x, double *location, int *nloc, double *scale, int *nscale, int *give_log)
{ /* pdf of univariate Laplace distribution */
  int nobs = *n, na = *nloc, nb = *nscale, log_pdf = *give_log;

  for (int i = 0; i < nobs; i++)
    y[i] = pdf_laplace(x[i], location[i % na], scale[i % nb], log_pdf);
}

void p_laplace(int *n, double *y, double *x, double *location, int *nloc, double *scale, int *nscale, int *lower_tail, int *log_p)
{ /* cdf of univariate Laplace distribution */
  int nobs = *n, na = *nloc, nb = *nscale, lower = *lower_tail, log_cdf = *log_p;

  for (int i = 0; i < nobs; i++)
    y[i] = cdf_laplace(x[i], location[i % na], scale[i % nb], lower, log_cdf);
}

void q_laplace(int *n, double *y, double *p, double *location, int *nloc, double *scale, int *nscale, int *lower_tail, int *log_p)
{ /* quantile function of univariate Laplace distribution */
  int nobs = *n, na = *nloc, nb = *nscale, lower = *lower_tail, log_prob = *log_p;

  for (int i = 0; i < nobs; i++)
    y[i] = quantile_laplace(p[i], location[i % na], scale[i % nb], lower, log_prob);
}

/* ========================================================================== *
 * pdf for the multivariate Laplace distribution
 * ========================================================================== */

void
pdf_mlaplace(double *y, double *x, int *nobs, int *vars, double *center, double *Scatter)
{ /* evaluation of the multivariate Laplace log-density function */
  int errcode = 0, job = 0, n = *nobs, p = *vars;
  double *Root, *z, D2, log_pdf;

  Root = (double *) Calloc(p * p, double);
  z    = (double *) Calloc(p, double);

  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in pdf_mlaplace gave code %d", errcode);
  
  log_pdf  = lgammafn(0.5 * p) - (double) p * M_LN_SQRT_PI - lgammafn((double) p) - (p + 1.0) * M_LN2;
  log_pdf -= logAbsDet(Root, p, p);
  
  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    D2 = mahalanobis(z, p, center, Root);
    y[i] = log_pdf - 0.5 * sqrt(D2);
  }

  Free(Root); Free(z);
}
