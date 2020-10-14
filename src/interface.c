/* ID: interface.c, last updated 2020-10-09, F.Osorio */

#include "base.h"
#include "interface.h"
#include <fastmatrix_API.h>

/* basic matrix manipulations */

void
ax_plus_y(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y */
  BLAS1_axpy(alpha, x, incx, y, incy, n);
}

void
scale(double *x, int n, int inc, double alpha)
{ /* x <- alpha * x */
  BLAS1_scale(alpha, x, inc, n);
}

double
sum_abs(double *x, int n, int inc)
{ /* sum(abs(x)) */
  return BLAS1_sum_abs(x, inc, n);
}

void
copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- x[,] */
  FM_copy_mat(y, ldy, x, ldx, nrow, ncol);
}

void
GAXPY(double *y, double alpha, double *a, int lda, int nrow, int ncol, double *x, double beta)
{ /* y <- alpha * a %*% x + beta * y */
  FM_GAXPY(y, alpha, a, lda, nrow, ncol, x, beta, 0); /* job = 0 */
}

void
mult_triangular_vec(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc)
{ /* wrapper to BLAS-2 routine 'DTRMV' */
  BLAS2_trmv(a, lda, n, uplo, trans, diag, x, inc);
}

void
mult_triangular_mat(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy)
{ /* wrapper to BLAS-3 routine 'DTRMM' */
  BLAS3_trmm(alpha, a, lda, nrow, ncol, side, uplo, trans, diag, y, ldy);
}

/* routines for matrix decompositions */

void
chol_decomp(double *a, int lda, int p, int job, int *info)
{ /* cholesky factorization of a real symmetric positive definite matrix a.
   * the factorization has the form:
   * a <- l %*% t(l) (job = 0), or a <- t(u) %*% u (job = 1),
   * where u is an upper triangular matrix and l is lower triangular */
  FM_chol_decomp(a, lda, p, job, info);
}

void
QR_decomp(double *mat, int ldmat, int nrow, int ncol, double *qraux, int *info)
{ /* return the QR decomposition of a rectangular matrix */
  FM_QR_decomp(mat, ldmat, nrow, ncol, qraux, info);
}

/* orthogonal-triangular operations */

void
QR_qty(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- qr.qty(qr, ymat) */
  FM_QR_qty(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, info);
}

void
QR_qy(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- qr.qy(qr, ymat) */
  FM_QR_qy(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, info);
}

/* triangular solver */

void
backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info)
{ /* backsolve solve triangular systems of the form r %*% x = b, where r
   * is an upper triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  FM_backsolve(r, ldr, n, b, ldb, nrhs, info);
}
