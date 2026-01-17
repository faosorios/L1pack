/* ID: R_init_L1pack.c, last updated 2025-06-23, F.Osorio */

#include "base.h"
#include "lad.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, nargs)  {#name, (DL_FUNC) &name, nargs}
#define F77DEF(name, nargs)   {#name, (DL_FUNC) &F77_NAME(name), nargs}

/* L1 estimation for linear regression */
extern void F77_NAME(l1)(int *, int *, int *, int *, double *, double *, double *, double *, double *, int *);
extern void lad_BR(double *, double *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);
extern void lad_EM(double *, double *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, int *);

/* multivariate Laplace estimation */
extern void Laplace_fitter(double *, int *, int *, double *, double *, double *, double *, double *, double *, int *);
extern void fitter_EQUAL(double *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *);

/* kernel U-statistics for asymptotic variance estimation of CCC */
void F77_NAME(rho1_ustat)(double *, double *, int *, double *, double *);

/* generalized spatial median */
extern void spatial_median(double *, int *, int *, double *, double *, double *, double *, double *, double *, int *, int *);

/* dpqr-functions for the Laplace distribution */
extern void d_laplace(int *, double *, double *, double *, int *, double *, int *, int *);
extern void p_laplace(int *, double *, double *, double *, int *, double *, int *, int *, int *);
extern void q_laplace(int *, double *, double *, double *, int *, double *, int *, int *, int *);
extern void r_laplace(int *, double *, double *, int *, double *, int *);

/* density and RNG for the multivariate Laplace distribution */
extern void pdf_mlaplace(double *, double *, int *, int *, double *, double *);
extern void RNG_mlaplace(double *, int *, double *, double *);

/* Wilson-Hilferty transformation */
extern void Wilson_Hilferty_Laplace(double *, int *, int *, double *);

/* registering C and symbols */
static const R_CMethodDef CEntries[]  = {
  CALLDEF(d_laplace,                8),
  CALLDEF(p_laplace,                9),
  CALLDEF(q_laplace,                9),
  CALLDEF(r_laplace,                6),
  CALLDEF(lad_BR,                  14),
  CALLDEF(lad_EM,                  13),
  CALLDEF(Laplace_fitter,          10),
  CALLDEF(fitter_EQUAL,            11),
  CALLDEF(pdf_mlaplace,             6),
  CALLDEF(RNG_mlaplace,             4),
  CALLDEF(spatial_median,          11),
  CALLDEF(Wilson_Hilferty_Laplace,  4),
  {NULL, NULL, 0}
};

/* registering F77 symbols */
static const R_FortranMethodDef F77Entries[] = {
  F77DEF(l1,              10),
  F77DEF(rho1_ustat,       5),
  {NULL, NULL, 0}
};

void R_init_L1pack(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, F77Entries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
