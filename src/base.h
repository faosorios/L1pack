/* ID: base.h, last updated 2024-05-17, F.Osorio */

#ifndef L1PACK_BASE_H
#define L1PACK_BASE_H

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rconfig.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

/* some definitions */
#define ABSTOL      1.0e-2
#define NULLP       (void *) 0
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))
#define MIN(a,b)    (((a)<(b)) ? (a) : (b))
#define SQR(x)      R_pow_di(x, 2)
#define DOUBLE_EPS  DBL_EPSILON
#define repeat      for(;;)

/* dims structure */
typedef struct DIMS_struct {
  int
    N,    /* total number of cases */
    n,    /* number of observations */
    p;    /* number of coefficients */
} DIMS_struct, *DIMS;

#endif /* L1PACK_BASE_H */
