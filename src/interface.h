/* ID: interface.h, last updated 2026-01-25, F.Osorio */

#ifndef L1PACK_INTERFACE_H
#define L1PACK_INTERFACE_H

/* basic matrix manipulations */
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void copy_mat(double *, int, double *, int, int, int);
extern void copy_lower(double *, int, double *, int, int);
extern void copy_vec(double *, int, double *, int, int);
extern double dot_product(double *, int, double *, int, int);
extern void GAXPY(double *, double, double *, int, int, int, double *, double);
extern double logAbsDet(double *, int, int);
extern void mult_triangular_vec(double *, int, int, char *, char *, char *, double *, int);
extern void mult_triangular_mat(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
extern double norm_sqr(double *, int, int);
extern double norm_two(double *, int, int);
extern void rank1_update(double *, int, int, int, double, double *, double *);
extern void scale(double *, int, int, double);
extern void solve_triangular_mat(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
extern void setzero(double *, int, int, int);
extern double sum_abs(double *, int, int);

/* routines for matrix decompositions */
extern void chol_decomp(double *, int, int, int, int *);
extern void QR_decomp(double *, int, int, int, double *, int *);

/* orthogonal-triangular operations */
extern void QR_qy(double *, int, int, int, double *, double *, int, int, int, int *);
extern void QR_qty(double *, int, int, int, double *, double *, int, int, int, int *);

/* triangular solver */
extern void backsolve(double *, int, int, double *, int, int, int *);

/* descriptive statistics */
extern void center_and_Scatter(double *, int, int, double *, double *, double *);
extern void center_online(double *, int, int, double *, double *);
extern void mediancenter(double *, int, int, double *, int *);
extern double sample_quantile(double *, int, int);

/* Mahalanobis distance */
extern double mahalanobis(double *, int, double *, double *);

/* Wilson-Hilferty transformation */
extern void WH_Laplace(double *, int, int, double *);

/* Misc */
extern void print_mat(double *, int, int, int, char *);

#endif /* L1PACK_INTERFACE_H */
