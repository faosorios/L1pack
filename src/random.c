#include "random.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);
/* ..end declarations */

/* 'dims' functions */

static DIMS
dims(int *pdims)
{ /* dims object */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->n = (int) pdims[0];
  ans->p = (int) pdims[1];
  return ans;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

/* multivariate Laplace random generation */

void
rand_laplace(double *y, int *pdims, double *center, double *Scatter)
{ /* multivariate Laplace random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int i, inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("DPOTRF in cholesky decomposition gave code %d", info);
  rand_spherical_laplace(y, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter,
    &(dm->p), y, &(dm->p));
  for (i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_laplace(double *y, int n, int p)
{ /* standard Laplace exponential variates */
  int i, j, one = 1;
  double radial;

  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++)
      y[j] = norm_rand();
    radial = sqrt(exp_rand());
    F77_CALL(dscal)(&p, &radial, y, &one);
    y += p;
  }
}
