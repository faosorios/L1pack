/* ID: lad.h, last updated 2023-12-26, F.Osorio */

#ifndef L1PACK_LAD_H
#define L1PACK_LAD_H

#include "base.h"

/* log-likelihood function for LAD regression */
double lad_logLik(double *, int);

/* routines for iterative reweighted least squares (IRLS)*/
int IRLS(double *, double *, int, int, double *, double *, double *, double *, double *, double *, int, double);
void IRLS_increment(double *, double *, int, int, double *, double *, double *, double *, double *, double *);

/* linear programming method to solve L1 regression problems */
void F77_NAME(l1br)(double *, double *, int *, int *, int *, int *, double *, double *, double *, int *, double *, int *, int *, int *);

/* standard Laplace random generation */
void rmlaplace_std(double *, int, int);

#endif /* L1PACK_LAD_H */
