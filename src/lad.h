#ifndef L1PACK_LAD_H
#define L1PACK_LAD_H

#include "base.h"
#include "matrix.h"

/* available methods */
typedef enum {
    BR,
    EM
} METHOD;

/* structure to hold model results */
typedef struct LAD_struct {
    DIMS dd;        /* dimension data info */
    double
        *y,         /* responses */
        *x,         /* model matrix */
        *coef,      /* coefficients estimates */
        *scale,     /* scale estimate */
        *sad,       /* minimum total absolute deviation */
        *logLik,    /* log-likelihood */
        *fitted,    /* fitted values */
        *resid,     /* residuals */
        *weights,   /* 'robust' weights computed for the EM algorithm */
        *settings,  /* additional settings */
        *control;   /* control settings for estimation algorithm */
    int
        maxIter;    /* maximun number of iterations */
    METHOD
        method;     /* algorithm of estimation */
    double
        tolerance;  /* convergence tolerance */
} LAD_struct, *LAD;

/* L1 estimation for linear models (to be called by R) */
extern void F77_NAME(l1)(int *, int *, int *, int *, double *, double *, double *, double *, double *, int *);
extern void lad(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void lad_acov(double *, int *, double *);

/* routines for estimation in linear models */
LAD lad_init(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *);
void lad_free(LAD);
void lad_fit(LAD);
double lad_objective(double *, int);
double lad_logLik(double *, int);

/* routines for Iterative Weighted Least Squares (IRLS)*/
int IRLS(double *, double *, DIMS, double *, double *, double *, double *, double *, double *, int, double);
void IRLS_increment(double *, double *, DIMS, double *, double *, double *, double *, double *, double *);
void qr_fitted(DIMS, double *, double *, double *, double *, double *);

/* wrapper for Fortran subroutine 'l1fit' */
int l1_fit(double *, double *, DIMS, double *, double *, double *, double *, double *, double *, double);
void F77_NAME(l1fit)(double *, double *, int *, int *, int *, int *, double *, double *, double *, int *, double *, int *, int *, int *);

#endif /* L1PACK_LAD_H */
