#ifndef L1PACK_MATRIX_H
#define L1PACK_MATRIX_H

#include "base.h"

/* basic matrix manipulations */
extern double sum_abs(double *, int, int);
extern double norm_sqr(double *, int, int);
extern double dot_product(double *, int, double *, int, int);
extern void scale(double *, int, int, double);
extern void zero_mat(double *, int, int, int);
extern void gaxpy(double *, double, double *, int, int, int, double *, double);
extern void copy_mat(double *, int, double *, int, int, int);
extern void add_mat(double *, int, double, double *, int, int, int);
extern void lower_tri(double *, int, double *, int, int, int);
extern void upper_tri(double *, int, double *, int, int, int);
extern void scale_mat(double *, int, double *, int, int, int, double);
extern void upper_mult_vec(double *, int, int, int, double *, double *);
extern void mult_mat(double *, int, int, int, double *, int, int, int, double *);
extern void crossprod(double *, int, int, int, double *, int, int, int, double *);
extern void outerprod(double *, int, int, int, double *, int, int, int, double *);
extern void rank1_update(double *, int, int, int, double, double *, double *);
extern double logAbsDet(double *, int, int);

/* routines for matrix decompositions */
extern void chol_decomp(double *, int, int, int, int *);
extern void svd_decomp(double *, int, int, int, double *, double *);
extern QRStruct QR_decomp(double *, int, int, int, double *);
extern void QR_free(QRStruct);

/* orthogonal-triangular operations */
extern void QR_qty(QRStruct, double *, int, int, int);
extern void QR_qy(QRStruct, double *, int, int);
extern void QR_store_R(QRStruct, double *, int);

/* matrix inversion and backsolve */
extern void invert_mat(double *, int, int);
extern void invert_triangular(int, double *, int, int);
extern void backsolve(int, double *, int, int, double *, int, int);

/* linear least-squares fit */
extern void lsfit(double *, int, int, int, double *, int, int, double *);

#endif /* L1PACK_MATRIX_H */
