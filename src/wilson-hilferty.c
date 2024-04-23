/* ID: wilson-hilferty.c, last updated 2024-01-03, F.Osorio */

#include "base.h"
#include "interface.h"

void
Wilson_Hilferty_Laplace(double *distances, int *n, int *p, double *z)
{ /* Wilson-Hilferty transformation for chi-squared variables */
  WH_Laplace(distances, *n, *p, z);
}
