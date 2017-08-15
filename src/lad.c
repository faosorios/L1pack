#include "lad.h"

/* static functions.. */
DIMS dims(int *);
void dims_free(DIMS);
double do_weight(double, double);
/* ..end declarations */

DIMS
dims(int *pdims)
{   /* dims object for linear models */
    DIMS ans;

    ans = (DIMS) Calloc(1, DIMS_struct);
    ans->n = (int) pdims[0];
    ans->p = (int) pdims[1];
    ans->N = ans->n * ans->p;
    return ans;
}

void
dims_free(DIMS this)
{   /* destructor for a dims object */
    Free(this);
}

double
do_weight(double residual, double sad)
{   /* compute weights */
    double ans, eps = R_pow(DOUBLE_EPS, .5);

    if (fabs(residual) < eps * sad)
        ans = 1.0;
    else
        ans = eps / fabs(residual);
    return ans;
}

void
lad(double *y, double *x, int *pdims, double *coef, double *scale, double *fitted,
    double *resid, double *weights, double *control, double *sad, double *logLik)
{   /* fitter for linear models under Laplace errors */
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
            iter = l1_fit(model->y, model->x, model->dd, model->coef, model->scale,
                          model->sad, model->fitted, model->resid, model->weights,
                          model->tolerance);
            break;
        case EM:
            iter = IRLS(model->y, model->x, model->dd, model->coef, model->scale,
                        model->sad, model->fitted, model->resid, model->weights,
                        model->maxIter, model->tolerance);
            break;
    }

    /* save numIter */
    (model->control)[3] = (double) iter;
}

LAD
lad_init(double *y, double *x, int *pdims, double *coef, double *scale, double *sad,
    double *fitted, double *resid, double *weights, double *control)
{   /* constructor for a linear model object */
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
{   /* destructor for a model object */
    dims_free(this->dd);
    Free(this);
}

double
lad_objective(double *residuals, int n)
{   /* total absolute deviation */
    int one = 1;
    double ans;

    ans = F77_CALL(dasum)(&n, residuals, &one);
    return ans;
}

double
lad_logLik(double *scale, int n)
{   /* log-likelihood function under Laplace errors */
    double ans;

    ans = (double) n * (0.5 * M_LN2 + 1.0 + log(*scale));
    return -ans;
}

int
IRLS(double *y, double *x, DIMS dd, double *coef, double *scale, double *sad,
    double *fitted, double *residuals, double *weights, int maxIter, double tolerance)
{   /* iteratively reweighted LS algorithm */
    int i, iter;
    double conv, SAD, newSAD, *incr, *working;

    /* initialization */
    incr    = (double *) Calloc(dd->p, double);
    working = (double *) Calloc(dd->n, double);
    SAD = lad_objective(residuals, dd->n);

    /* main loop */
    for (iter = 1; iter <= maxIter; iter++) {
        /* E-step */
        for (i = 0; i < dd->n; i++)
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
{   /* increment for direction search in IRLS */
    int i, j, info = 0, one = 1;
    char *uplo = "U", *diag = "N", *side = "L", *notrans = "N", *trans = "T";
    double stepsize = 1.0, wts, *z, *qraux, *work;

    /* initialization */
    z     = (double *) Calloc(dd->n * dd->p, double);
    qraux = (double *) Calloc(dd->p, double);
    work  = (double *) Calloc(dd->p, double);

    /* transformed model matrix and working residuals */
    for (i = 0; i < dd->n; i++) {
        wts = sqrt(weights[i]);
        working[i] = wts * residuals[i];
        for (j = 0; j < dd->p; j++)
            z[i + j * dd->n] = wts * x[i + j * dd->n];
    }

    /* solve the transformed LS-problem */
    F77_CALL(dgeqrf)(&(dd->n), &(dd->p), z, &(dd->n), qraux, work, &(dd->p),
                     &info);
    if (info)
        error("DGEQRF in IRLS_increment gave code %d", info);
    F77_CALL(dormqr)(side, trans, &(dd->n), &one, &(dd->p), z, &(dd->n), qraux,
                     working, &(dd->n), work, &(dd->p), &info);
    if (info)
        error("DORMQR in IRLS_increment gave code %d", info);
    Memcpy(incr, working, dd->p);
    F77_CALL(dtrtrs)(uplo, notrans, diag, &(dd->p), &one, z, &(dd->n), incr,
                     &(dd->p), &info);
    if (info)
        error("DTRTRS in IRLS_increment gave code %d", info);
    /* update coefficients */
    F77_CALL(daxpy)(&(dd->p), &stepsize, incr, &one, coef, &one);

    /* fitted values and residuals */
    qr_fitted(dd, z, coef, fitted, qraux, work);
    for (i = 0; i < dd->n; i++) {
        wts = sqrt(weights[i]);
        fitted[i] /= wts;
        residuals[i] = y[i] - fitted[i];
    }
    Free(z); Free(qraux); Free(work);
}

void
qr_fitted(DIMS dd, double *z, double *coef, double *fitted, double *qraux,
    double *work)
{   /* compute the fitted values */
    int i, info = 0, one = 1;
    char *uplo = "U", *diag = "N", *side = "L", *notrans = "N";

    for (i = 0; i < dd->n; i++)
        fitted[i] = 0.0;
    Memcpy(fitted, coef, dd->p);
    F77_CALL(dtrmv)(uplo, notrans, diag, &(dd->p), z, &(dd->n), fitted, &one);
    F77_CALL(dormqr)(side, notrans, &(dd->n), &one, &(dd->p), z, &(dd->n),
                     qraux, fitted, &(dd->n), work, &(dd->p), &info);
    if (info)
        error("DORMQR in qr_fitted gave code %d", info);
}

int
l1_fit(double *y, double *x, DIMS dd, double *coef, double *scale, double *sad,
    double *fitted, double *residuals, double *weights, double tolerance)
{   /* wrapper for 'l1fit' */
    int i, n2 = dd->n + 2, p2 = dd->p + 2, info, iter, rank, *work;
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
    gaxpy(fitted, 1.0, x, dd->n, dd->n, dd->p, coef, 1.0);
    for (i = 0; i < dd->n; i++) {
        /* flag for 'basic' observations */
        if (residuals[i] == 0.0)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    if (info == 0)
        error("L1FIT optimal solution is probably non-unique.");
    if (info != 1)
        error("L1FIT calculations terminated prematurely");
    if (rank != dd->p)
        error("L1FIT matrix not of full rank, apparent rank %d", rank);

    return iter;
}

void
lad_acov(double *R, int *pdims, double *acov)
{ /* unscaled asymptotic covariance matrix */
  int info = 0, job = 1;
  DIMS dd;

  /* unscaled asymptotic covariance */
  dd = dims(pdims);
  invert_triangular(job, R, dd->p, dd->p);
  if (info)
    error("DTRTRI in lad_acov gave code %d", info);
  outerprod(R, dd->p, dd->p, dd->p, R, dd->p, dd->p, dd->p, acov);

  dims_free(dd);
}
