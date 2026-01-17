/* ID: fitter_equal.c, last updated 2025-06-18, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

/* static functions.. */
static double update_lambda(int, double *, double *);
static void update_Scatter(double *, int, int, double *, double *, double *);
/* ..end declarations */

void
fitter_EQUAL(double *x, int *nobs, int *vars, double *center, double *lambda, double *Scatter, 
  double *distances, double *weights, double *logLik, double *tolerance, int *maxiter)
{ /* restricted estimation for the multivariate Laplace distribution under homogeneity 
   * of means, described in Appendix of Vallejos, Osorio & Ferrer (2025) */
  int errcode = 0, iter = 0, job = 0, n = *nobs, p = *vars, maxit = *maxiter, npars;
  double conv, fnc = *logLik, newfnc, *mean, *Root, tol = *tolerance;

  mean = (double *) R_Calloc(p, double);
  Root = (double *) R_Calloc(p * p, double);

  /* Cholesky decomposition of Scatter matrix */
  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in Laplace fitter gave code %d", errcode);

  /* set elements of the restricted mean */
  for (int j = 0; j < p; j++)
    mean[j] = *lambda;
  
  /* main loop */
  repeat {
    /* E-step */
    E_step(x, n, p, mean, Root, distances, weights);

    /* M-step */
    center_online(x, n, p, weights, center);
    *lambda = update_lambda(p, center, Root);
    update_Scatter(x, n, p, weights, lambda, Scatter);

    /* 'update' mean */
    for (int j = 0; j < p; j++)
      mean[j] = *lambda;

    iter++;

    /* evaluating objective function */
    copy_lower(Root, p, Scatter, p, p);
    chol_decomp(Root, p, p, job, &errcode);
    if (errcode)
      error("Cholesky decomposition in Laplace fitter gave code %d", errcode);
    newfnc = logLik_Laplace(distances, n, p, Root);
    
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
    error("Cholesky decomposition in Laplace fitter gave code %d", errcode);
  *logLik = logLik_Laplace(distances, n, p, Root);

  R_Free(mean); R_Free(Root);
}

static double
update_lambda(int p, double *center, double *Root)
{ /* update lambda (common factor of the location) estimate */
  char *uplo = "L", *notrans = "N", *diag = "N";
  double lambda, *ones, *z;

  ones = (double *) R_Calloc(p, double);
  z    = (double *) R_Calloc(p, double);

  for (int j = 0; j < p; j++) {
    ones[j] = 1.0;
    z[j] = center[j];
  }

  mult_triangular_vec(Root, p, p, uplo, notrans, diag, ones, 1);
  mult_triangular_vec(Root, p, p, uplo, notrans, diag, z, 1);

  lambda  = dot_product(ones, 1, z, 1, p);
  lambda /= dot_product(ones, 1, ones, 1, p);

  R_Free(ones); R_Free(z);

  return lambda;
}

static void
update_Scatter(double *x, int n, int p, double *weights, double *lambda, double *Scatter)
{ /* compute the restricted Scatter estimate */
  double *z, wts;

  /* initialization */
  z = (double *) R_Calloc(p, double);
  setzero(Scatter, p, p, p);

  /* updating stage */
  for (int i = 0; i < n; i++) {
    wts = weights[i];
    copy_vec(z, 1, x + i, n, p);
    for (int j = 0; j < p; j++)
      z[j] -= *lambda;
    rank1_update(Scatter, p, p, p, wts / n, z, z);
  }

  R_Free(z);
}
