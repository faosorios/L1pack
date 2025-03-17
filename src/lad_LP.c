/* ID: lad_LP.c, last updated 2024-09-07, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

void
lad_BR(double *y, double *x, int *nobs, int *vars, double *coef, double *scale,
  double *sad, double *fitted, double *residuals, double *logLik, double *tolerance,
  int *rank, int *iter, int *info)
{ /* lad_BR uses a modification of the simplex method of linear programming to
   * calculate an L1 solution to an over-determined system of linear equations.
   * Algorithm 478: Commun. ACM 17, 1974, 319-320. doi: 10.1145/355616.361024
   * This is a wrapper for Fortran 'L1BR' routine */
  int n = *nobs, p = *vars, n2, p2, *work;
  double minimum, *z;

  /* initialization */
  n2   = n + 2;
  p2   = p + 2;
  z    = (double *) R_Calloc(n2 * p2, double);
  work = (int *) R_Calloc(n, int);
  copy_mat(z, n2, x, n, n, p);

  /* call fitter */
  F77_CALL(l1br)(z, y, &n, &p, &n2, &p2, coef, residuals, &minimum,
                 iter, tolerance, rank, info, work);
  R_Free(z); R_Free(work);

  /* save results */
  *sad   = minimum;
  *scale = M_SQRT2 * minimum / n;

  /* compute fitted values */
  GAXPY(fitted, 1.0, x, n, n, p, coef, 1.0);

  /* log-likelihood evaluation */
  *logLik = lad_logLik(scale, n);
}
