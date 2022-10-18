/* ID: lad.h, last updated 2022-10-17, F.Osorio */

#ifndef L1PACK_LAD_H
#define L1PACK_LAD_H

#include "base.h"

/* routines for LAD estimation in linear regression */
double do_weight(double, double);
double lad_objective(double *, int);
double lad_logLik(double *, int);

/* routines for iterative reweighted least squares (IRLS)*/
int IRLS(double *, double *, int, int, double *, double *, double *, double *, double *, double *, int, double);
void IRLS_increment(double *, double *, int, int, double *, double *, double *, double *, double *, double *);

/* linear programming method to solve L1 regression problems */
void F77_NAME(l1br)(double *, double *, int *, int *, int *, int *, double *, double *, double *, int *, double *, int *, int *, int *);

/* spherical random generation */
void rand_spherical_laplace(double *, int, int);

#endif /* L1PACK_LAD_H */
