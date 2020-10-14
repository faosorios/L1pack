/* ID: random.c, last updated 2020-10-08, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

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
  int info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();

  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("DPOTRF in rand_laplace gave error code %d", info);

  rand_spherical_laplace(y, dm->n, dm->p);
  mult_triangular_mat(1.0, Scatter, dm->p, dm->p, dm->n, side, uplo, trans, diag, y, dm->p);

  for (int i = 0; i < dm->n; i++) {
    ax_plus_y(1.0, center, 1, y, 1, dm->p);
    y += dm->p;
  }

  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_laplace(double *y, int n, int p)
{ /* standard Laplace exponential variates */
  double radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    radial = sqrt(exp_rand());
    scale(y, p, 1, radial);
    y += p;
  }
}
