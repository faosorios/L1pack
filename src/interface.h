/* ID: interface.h, last updated 2020-10-09, F.Osorio */

#ifndef L1PACK_INTERFACE_H
#define L1PACK_INTERFACE_H

/* basic matrix manipulations */
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void copy_mat(double *, int, double *, int, int, int);
extern void GAXPY(double *, double, double *, int, int, int, double *, double);
extern void mult_triangular_vec(double *, int, int, char *, char *, char *, double *, int);
extern void mult_triangular_mat(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
extern void scale(double *, int, int, double);
extern double sum_abs(double *, int, int);

/* routines for matrix decompositions */
extern void chol_decomp(double *, int, int, int, int *);
extern void QR_decomp(double *, int, int, int, double *, int *);

/* orthogonal-triangular operations */
extern void QR_qy(double *, int, int, int, double *, double *, int, int, int, int *);
extern void QR_qty(double *, int, int, int, double *, double *, int, int, int, int *);

/* triangular solver */
extern void backsolve(double *, int, int, double *, int, int, int *);

#endif /* L1PACK_INTERFACE_H */
