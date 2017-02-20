#ifndef L1PACK_BASE_H
#define L1PACK_BASE_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>

/* some definitions */
#define NULLP    (void *) 0
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SQR(x)   R_pow_di(x, 2)
#define ETA_CONV 1.0e-2
#define repeat for(;;)

/* dims structure */
typedef struct DIMS_struct {
    int
        N,      /* total number of cases */
        n,      /* number of observations */
        p;      /* number of coefficients */
} DIMS_struct, *DIMS;

/* QR structure */
typedef struct QR_struct {
    double *mat, *qraux;
    int ldmat, nrow, ncol;
} QR_struct, *QRStruct;

#endif /* L1PACK_BASE_H */
