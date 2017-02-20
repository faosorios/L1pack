#ifndef L1PACK_RANDOM_H
#define L1PACK_RANDOM_H

#include "matrix.h"

/* multivariate Laplace random generation (to be called by R) */
extern void rand_laplace(double *, int *, double *, double *);

/* spherical random generation */
extern void rand_spherical_laplace(double *, int, int);

#endif /* L1PACK_RANDOM_H */
