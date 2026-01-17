/* ID: laplace_fitter.c, last updated 2024-09-07, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

double
do_weight(double length, double distance) 
{ /* multivariate Laplace weight */
  double d, ratio, wts, x;

  d = sqrt(distance);
  x = 0.5 * d;
  ratio = bessel_k(x, 0.5 * length - 1.0, 1.0) / bessel_k(x, 0.5 * length, 1.0);
  wts = 0.5 * ratio / d;

  return wts;
}

void
Laplace_fitter(double *x, int *nobs, int *vars, double *center, double *Scatter, double *distances, 
  double *weights, double *logLik, double *tolerance, int *maxiter)
{ /* fits the multivariate Laplace model considering an unstructured covariance matrix 
   * based on the EM algorithm proposed by Yavuz & Arslan (2018). Stat Papers 59, 271-289. 
   * doi: 10.1007/s00362-016-0763-x */
  int errcode = 0, iter = 0, job = 0, n = *nobs, p = *vars, maxit = *maxiter;
  double conv, fnc = *logLik, newfnc, *Root, tol = *tolerance;

  Root  = (double *) R_Calloc(p * p, double);

  /* Cholesky decomposition of Scatter matrix */
  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in Laplace fitter gave code %d", errcode);
  
  /* main loop */
  repeat {
    /* E-step */
    E_step(x, n, p, center, Root, distances, weights);

    /* M-step */
    center_and_Scatter(x, n, p, weights, center, Scatter);

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

  R_Free(Root);
}

void
E_step(double *x, int n, int p, double *center, double *Root, double *distances, double *weights)
{ /* computation of Mahalanobis distances and weights for the Laplace distribution */
  double *z;

  z = (double *) R_Calloc(p, double);

  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    distances[i] = mahalanobis(z, p, center, Root);
    weights[i] = do_weight((double) p, distances[i]);
  }

  R_Free(z);
}

double 
logLik_Laplace(double *distances, int n, int p, double *Root)
{ /* evaluate the Laplace log-likelihood */
  double accum = 0.0, val;

  /* sum the kernel of the log-density */
  for (int i = 0; i < n; i++)
    accum += sqrt(*distances++);

  val  = lgammafn(0.5 * p) - (double) p * M_LN_SQRT_PI - lgammafn((double) p) - (p + 1.0) * M_LN2;
  val -= logAbsDet(Root, p, p);
  val *= (double) n;
  val -= 0.5 * accum;

  return val;
}
