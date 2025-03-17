/* ID: spatial_median.c, last updated 2024-09-07, F.Osorio */

#include "base.h"
#include "interface.h"

/* static functions.. */
static double do_weight(double);
static void stage_center(double *, int, int, double *, double *, int *);
static void stage_weights(double *, int, int, double *, double *, double *, double *);
static void stage_Scatter(double *, int, int, double *, double *, double *);
static double logLik_Kotz(double *, int, int, double *);
/* ..end declarations */

static double
do_weight(double distance) 
{ /* Kotz-type weight */
  double wts;
  wts = 1.0 / sqrt(distance);
  return wts;
}

void
spatial_median(double *x, int *nobs, int *vars, double *median, double *Scatter, double *distances, 
  double *weights, double *logLik, double *tolerance, int *maxiter, int *iterations)
{ /* fits the multivariate Laplace model considering an unstructured covariance matrix */
  int errcode = 0, iter = 0, inner = 0, job = 0, n = *nobs, p = *vars, maxit = *maxiter;
  double conv, fnc = *logLik, newfnc, *Root, tol = *tolerance;

  Root = (double *) R_Calloc(p * p, double);

  /* Cholesky decomposition of Scatter matrix */
  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in spatial_median gave code %d", errcode);
  
  /* main loop */
  repeat {
    /* median update */
    stage_center(x, n, p, median, Root, &inner);
    *iterations += inner;

    /* Scatter update */
    stage_weights(x, n, p, median, Root, distances, weights);
    stage_Scatter(x, n, p, weights, median, Scatter);

    iter++;

    /* evaluating objective function */
    copy_lower(Root, p, Scatter, p, p);
    chol_decomp(Root, p, p, job, &errcode);
    if (errcode)
      error("Cholesky decomposition in spatial_median gave code %d", errcode);
    newfnc = logLik_Kotz(distances, n, p, Root);
    
    /* eval convergence */
    conv = fabs((newfnc - fnc) / (newfnc + ABSTOL));
    if (conv < tol)
      break; /* successful completion */
    if (iter >= maxit)
      break;  /* maximum number of iterations exceeded */

    fnc = newfnc;
  }

  *maxiter = iter;

  /* computation of log-likelihood function */
  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in spatial_median gave code %d", errcode);
  *logLik = logLik_Kotz(distances, n, p, Root);

  R_Free(Root);
}

static void
stage_center(double *x, int n, int p, double *median, double *Root, int *iterations)
{ /* computation of generalized spatial median (or median center) */
  char *side = "R", *lower = "L", *trans = "T", *notrans = "N", *diag = "N";
  double *z;

  z = (double *) R_Calloc(n * p, double);

  copy_mat(z, n, x, n, n, p);
  solve_triangular_mat(1.0, Root, p, n, p, side, lower, trans, diag, z, n);
  mediancenter(z, n, p, median, iterations);
  mult_triangular_vec(Root, p, p, lower, notrans, diag, median, 1);
  R_Free(z);
}

static void
stage_weights(double *x, int n, int p, double *median, double *Root, double *distances, double *weights)
{ /* computation of Mahalanobis distances and weights for the Kotz-type distribution */
  double *z;

  z = (double *) R_Calloc(p, double);

  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    distances[i] = mahalanobis(z, p, median, Root);
    weights[i] = do_weight(distances[i]);
  }

  R_Free(z);
}

static void
stage_Scatter(double *x, int n, int p, double *weights, double *median, double *Scatter)
{ /* compute the restricted Scatter estimate */
  double *z, wts;

  /* initialization */
  z = (double *) R_Calloc(p, double);
  setzero(Scatter, p, p, p);

  /* updating stage */
  for (int i = 0; i < n; i++) {
    wts = weights[i];
    copy_vec(z, 1, x + i, n, p);
    ax_plus_y(-1.0, median, 1, z, 1, p);
    rank1_update(Scatter, p, p, p, wts / n, z, z);
  }

  R_Free(z);
}

static double 
logLik_Kotz(double *distances, int n, int p, double *Root)
{ /* evaluate the Kotz-type log-likelihood */
  double accum = 0.0, val;

  /* sum the kernel of the log-density */
  for (int i = 0; i < n; i++)
    accum += sqrt(*distances++);

  val  = -1. * p * M_LN2 - (p - 1.) * M_LN_SQRT_PI - lgammafn(0.5 * (p + 1.));
  val -= logAbsDet(Root, p, p);
  val *= (double) n;
  val -= accum;

  return val;
}
