/* ID: RNG.c, last updated 2023-05-17, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);
static double rand_laplace(double, double);
/* ..end declarations */

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

/* ========================================================================== *
 * univariate Laplace random generation
 * ========================================================================== */

void r_laplace(int *n, double *x, double *location, int *nloc, double *scale, int *nscale)
{ /* univariate Laplace random generation */
  int i, nobs = *n, na = *nloc, nb = *nscale;

  GetRNGstate();
  for (i = 0; i < nobs; i++)
    x[i] = rand_laplace(location[i % na], scale[i % nb]);
  PutRNGstate();
}

static double rand_laplace(double location, double scale)
{ /* random variates from the univariate Laplace distribution */
  double u;

  u = unif_rand();

  if (scale == 0.)
    return location;
  
  if (u >= .5)
    return location - M_SQRT1_2 * scale * log(2. * (1. - u));
  else
    return location + M_SQRT1_2 * scale * log(2. * u);
}

/* ========================================================================== *
 * multivariate Laplace random generation
 * ========================================================================== */

void
RNG_mlaplace(double *y, int *pdims, double *center, double *Scatter)
{ /* multivariate Laplace random number generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  int info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();

  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("DPOTRF in rand_laplace gave error code %d", info);

  rmlaplace_std(y, dm->n, dm->p);
  mult_triangular_mat(1.0, Scatter, dm->p, dm->p, dm->n, side, uplo, trans, diag, y, dm->p);

  for (int i = 0; i < dm->n; i++) {
    ax_plus_y(1.0, center, 1, y, 1, dm->p);
    y += dm->p;
  }

  PutRNGstate();
  dims_free(dm);
}

void
rmlaplace_std(double *y, int n, int p)
{ /* standard multivariate Laplace deviates */
  double radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    radial = rgamma((double) p, 2.);
    radial /= norm_two(y, p, 1);
    scale(y, p, 1, radial);
    y += p;
  }
}
