/* ID: lad_EM.c, last updated 2024-09-07, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

/* static functions.. */
static double weight_1D(double, double);
static double lad_objective(double *, int);
/* ..end declarations */

static double
weight_1D(double residual, double eps)
{ /* weighting strategy based on Phillips (2002) */
  double ans, dev;

  dev = fabs(residual);
  if (dev < eps)
    ans = 1.0;        /* basic observation */
  else
    ans = 1.0 / dev;  /* 'EM' weights */
  return ans;
}

void
lad_EM(double *y, double *x, int *nobs, int *vars, double *coef, double *scale,
  double *fitted, double *resid, double *weights, double *sad, double *logLik,
  double *tolerance, int *maxiter)
{ /* lad_EM uses an EM algorithm based on iteratively reweighted least
   * squares (IRLS) scheme to approximate the LAD solution.
   * Phillips (2002). Stat. Comput. 12, 281-285. doi: 10.1023/A:1020759012226 */
  int iter = 0, n = *nobs, p = *vars, maxit = *maxiter;
  double tol = *tolerance;

  iter = IRLS(x, y, n, p, coef, scale, sad, fitted, resid, weights, maxit, tol);
  *maxiter = iter;
  *logLik  = lad_logLik(scale, n);
}

static double
lad_objective(double *residuals, int n)
{ /* sum of absolute deviations */
  double ans;

  ans = sum_abs(residuals, n, 1);
  return ans;
}

double
lad_logLik(double *scale, int n)
{ /* log-likelihood function under Laplace errors */
  double ans;

  ans = (double) n * (0.5 * M_LN2 + 1.0 + log(*scale));
  return -ans;
}

int
IRLS(double *x, double *y, int n, int p, double *coef, double *scale, double *sad,
  double *fitted, double *residuals, double *weights, int maxiter, double tolerance)
{ /* iteratively reweighted LS algorithm */
  int iter;
  double conv, eps = R_pow(DOUBLE_EPS, .5), SAD, newSAD, *incr, *working;

  /* initialization */
  incr    = (double *) R_Calloc(p, double);
  working = (double *) R_Calloc(n, double);
  SAD = lad_objective(residuals, n);

  /* main loop */
  for (iter = 1; iter <= maxiter; iter++) {
    /* E-step */
    for (int i = 0; i < n; i++)
      weights[i] = weight_1D(residuals[i], eps);

    /* M-step */
    IRLS_increment(x, y, n, p, coef, incr, working, fitted, residuals, weights);
    newSAD = lad_objective(residuals, n);
    *sad = newSAD;
    *scale = M_SQRT2 * *sad / n; /* CM-step for the 'scale' parameter */

    /* eval convergence */
    conv = fabs((newSAD - SAD) / (newSAD + ABSTOL));
    if (conv < tolerance) { /* successful completion */
      R_Free(incr); R_Free(working);
      return iter;
    }
    SAD = newSAD;
  }
  R_Free(incr); R_Free(working);

  return (iter - 1);
}

void
IRLS_increment(double *x, double *y, int n, int p, double *coef, double *incr,
  double *working, double *fitted, double *residuals, double *weights)
{ /* increment for direction search in IRLS */
  int info = 0;
  double stepsize = 1.0, wts, *z, *qraux;

  /* initialization */
  z     = (double *) R_Calloc(n * p, double);
  qraux = (double *) R_Calloc(p, double);

  /* transformed model matrix and working residuals */
  for (int i = 0; i < n; i++) {
    wts = sqrt(weights[i]);
    working[i] = wts * residuals[i];
    for (int j = 0; j < p; j++)
      z[i + j * n] = wts * x[i + j * n];
  }

  /* QR decomposition of transformed model matrix */
  QR_decomp(z, n, n, p, qraux, &info);
  if (info)
    error("QR_decomp in IRLS_increment gave error code %d", info);

  /* solve the transformed LS problem */
  QR_qty(z, n, n, p, qraux, working, n, n, 1, &info);
  if (info)
    error("QR_qty in IRLS_increment gave error code %d", info);
  Memcpy(incr, working, p);
  backsolve(z, n, p, incr, p, 1, &info);
  if (info)
    error("backsolve in IRLS_increment gave error code %d", info);

  /* update coefficients */
  ax_plus_y(stepsize, incr, 1, coef, 1, p);

  /* fitted values */
  for (int i = 0; i < n; i++)
    fitted[i] = 0.0;
  Memcpy(fitted, coef, p);
  mult_triangular_vec(z, n, p, "U", "N", "N", fitted, 1);
  QR_qy(z, n, n, p, qraux, fitted, n, n, 1, &info);
  if (info)
    error("QR_qy in IRLS_increment gave error code %d", info);

  /* fitted values and residuals in original scale */
  for (int i = 0; i < n; i++) {
    wts = sqrt(weights[i]);
    fitted[i] /= wts;
    residuals[i] = y[i] - fitted[i];
  }
  R_Free(z); R_Free(qraux);
}
