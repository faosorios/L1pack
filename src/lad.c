/* ID: lad.c, last updated 2020-10-09, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

/* static functions.. */
DIMS dims(int *);
void dims_free(DIMS);
double do_weight(double, double);
/* ..end declarations */

DIMS
dims(int *pdims)
{ /* dims object for linear models */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->n = (int) pdims[0];
  ans->p = (int) pdims[1];
  ans->N = ans->n * ans->p;
  return ans;
}

void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

double
do_weight(double residual, double sad)
{ /* compute weights */
  double ans, eps = R_pow(DOUBLE_EPS, .5);

  if (fabs(residual) < eps * sad)
    ans = 1.0;
  else
    ans = eps / fabs(residual);
  return ans;
}

void
lad_fitter(double *y, double *x, int *pdims, double *coef, double *scale, double *fitted,
  double *resid, double *weights, double *control, double *sad, double *logLik)
{ /* fitter for linear models under Laplace errors */
  LAD model;

  model = lad_init(y, x, pdims, coef, scale, sad, fitted, resid, weights, control);
  lad_fit(model);
  *logLik = lad_logLik(model->scale, model->dd->n);
  lad_free(model);
}

void
lad_fit(LAD model)
{
  int iter = 0;

  switch(model->method) {
    case BR:
      iter = l1fit_BR(model->y, model->x, model->dd, model->coef, model->scale, model->sad,
                      model->fitted, model->resid, model->weights, model->tolerance);
      break;
    case EM:
      iter = IRLS(model->y, model->x, model->dd, model->coef, model->scale, model->sad,
                  model->fitted, model->resid, model->weights, model->maxIter, model->tolerance);
      break;
  }

  /* save numIter */
  (model->control)[3] = (double) iter;
}

LAD
lad_init(double *y, double *x, int *pdims, double *coef, double *scale, double *sad,
  double *fitted, double *resid, double *weights, double *control)
{ /* constructor for a linear model object */
  LAD model;

  model = (LAD) Calloc(1, LAD_struct);
  model->dd = dims(pdims);
  model->y = y;
  model->x = x;
  model->coef = coef;
  model->scale = scale;
  model->sad = sad;
  model->fitted = fitted;
  model->resid = resid;
  model->weights = weights;
  model->control = control;
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->method = (int) control[2];
  return model;
}

void
lad_free(LAD this)
{ /* destructor for a model object */
  dims_free(this->dd);
  Free(this);
}

double
lad_objective(double *residuals, int n)
{ /* total absolute deviation */
  double ans;

  ans = sum_abs(residuals, n, 1);
  return ans;
}

double
lad_logLik(double *scale, int n)
{ /* log-likelihood function under Laplace errors */
  double ans;

  ans = (double) n * (0.5 * M_LN2 + 1.0 + log(*scale));
  return -ans;
}

int
IRLS(double *y, double *x, DIMS dd, double *coef, double *scale, double *sad, double *fitted,
  double *residuals, double *weights, int maxIter, double tolerance)
{ /* iteratively reweighted LS algorithm */
  int iter;
  double conv, SAD, newSAD, *incr, *working;

  /* initialization */
  incr    = (double *) Calloc(dd->p, double);
  working = (double *) Calloc(dd->n, double);
  SAD = lad_objective(residuals, dd->n);

  /* main loop */
  for (iter = 1; iter <= maxIter; iter++) {
    /* E-step */
    for (int i = 0; i < dd->n; i++)
      weights[i] = do_weight(residuals[i], SAD);

    /* M-step */
    IRLS_increment(y, x, dd, fitted, residuals, weights, coef, incr, working);
    newSAD = lad_objective(residuals, dd->n);
    *sad = newSAD;
    *scale = M_SQRT2 * *sad / dd->n;

    /* eval convergence */
    conv = fabs((newSAD - SAD) / (newSAD + ETA_CONV));
    if (conv < tolerance) { /* successful completion */
      Free(incr); Free(working);
      return iter;
    }
    SAD = newSAD;
  }
  Free(incr); Free(working);

  return (iter - 1);
}

void
IRLS_increment(double *y, double *x, DIMS dd, double *fitted, double *residuals,
  double *weights, double *coef, double *incr, double *working)
{ /* increment for direction search in IRLS */
  int info = 0;
  double stepsize = 1.0, wts, *z, *qraux;

  /* initialization */
  z     = (double *) Calloc(dd->n * dd->p, double);
  qraux = (double *) Calloc(dd->p, double);

  /* transformed model matrix and working residuals */
  for (int i = 0; i < dd->n; i++) {
    wts = sqrt(weights[i]);
    working[i] = wts * residuals[i];
    for (int j = 0; j < dd->p; j++)
      z[i + j * dd->n] = wts * x[i + j * dd->n];
  }

  /* QR decomposition of transformed model matrix */
  QR_decomp(z, dd->n, dd->n, dd->p, qraux, &info);
  if (info)
    error("QR_decomp in IRLS_increment gave error code %d", info);

  /* solve the transformed LS problem */
  QR_qty(z, dd->n, dd->n, dd->p, qraux, working, dd->n, dd->n, 1, &info);
  if (info)
    error("QR_qty in IRLS_increment gave error code %d", info);
  Memcpy(incr, working, dd->p);
  backsolve(z, dd->n, dd->p, incr, dd->p, 1, &info);
  if (info)
    error("DTRTRS in IRLS_increment gave error code %d", info);

  /* update coefficients */
  ax_plus_y(stepsize, incr, 1, coef, 1, dd->p);

  /* fitted values */
  for (int i = 0; i < dd->n; i++)
    fitted[i] = 0.0;
  Memcpy(fitted, coef, dd->p);
  mult_triangular_vec(z, dd->n, dd->p, "U", "N", "N", fitted, 1);
  QR_qy(z, dd->n, dd->n, dd->p, qraux, fitted, dd->n, dd->n, 1, &info);
  if (info)
    error("QR_qy in IRLS_increment gave error code %d", info);

  /* fitted values and residuals in original scale */
  for (int i = 0; i < dd->n; i++) {
    wts = sqrt(weights[i]);
    fitted[i] /= wts;
    residuals[i] = y[i] - fitted[i];
  }
  Free(z); Free(qraux);
}

int
l1fit_BR(double *y, double *x, DIMS dd, double *coef, double *scale, double *sad,
  double *fitted, double *residuals, double *weights, double tolerance)
{ /* wrapper for 'l1fit' */
  int n2 = dd->n + 2, p2 = dd->p + 2, info, iter, rank, *work;
  double minimum, *z;

  /* initialization */
  z    = (double *) Calloc(n2 * p2, double);
  work = (int *) Calloc(dd->n, int);
  copy_mat(z, n2, x, dd->n, dd->n, dd->p);

  /* call fitter */
  F77_CALL(l1fit)(z, y, &(dd->n), &(dd->p), &n2, &p2, coef, residuals, &minimum,
                  &iter, &tolerance, &rank, &info, work);
  Free(z); Free(work);

  /* save results */
  *sad   = minimum;
  *scale = M_SQRT2 * minimum / dd->n;

  /* post-processing */
  GAXPY(fitted, 1.0, x, dd->n, dd->n, dd->p, coef, 1.0);
  for (int i = 0; i < dd->n; i++) {
    /* flag for 'basic' observations */
    if (residuals[i] == 0.0)
      weights[i] = 1.0;
    else
      weights[i] = 0.0;
  }

  if (info == 0)
    error("l1fit: optimal solution is probably non-unique.");
  if (info != 1)
    error("l1fit: calculations terminated prematurely");
  if (rank != dd->p)
    error("l1fit: matrix not of full rank, apparent rank %d", rank);

  return iter;
}
