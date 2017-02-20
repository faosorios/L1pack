#ifndef L1PACK_DISTN_H
#define L1PACK_DISTN_H

#include "base.h"

/* dpqr-functions for Laplace distribution (to be called by R) */
extern void dlaplace(int *, double *, double *, double *, int *, double *, int *, int *);
extern void plaplace(int *, double *, double *, double *, int *, double *, int *, int *, int *);
extern void qlaplace(int *, double *, double *, double *, int *, double *, int *, int *, int *);
extern void rlaplace(int *, double *, double *, int *, double *, int *);

#endif /* L1PACK_DISTN_H */
