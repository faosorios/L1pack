/* ID: lad_nuisance.c, last updated 2026-01-25, F.Osorio */

#include "base.h"
#include "interface.h"
#include "lad.h"

void
nuisance_vcov(double *resid, int *nok, double *alpha, int *job, double *value)
{ /* estimation of the 'nuisance' parameter representing the asymptotic 
   * variance of the sample median from the distribution of the random 
   * disturbances */
  int lo, hi, m = *nok, task = *job;
  double ans, conf = *alpha, k, z;

  z = qnorm(1. - conf, 0.0, 1.0, 1, 0);
  ans = .5 * sqrt(m) / z;
  if (task) { 
    /* McKean and Schrader (1987), see description at
     * Gonin and Money (1989), pag. 19. doi: 10.1201/9780203745526 */
    k = (m + 1) / 2. - z * sqrt(m / 4.);
    lo = (int) k + .5;
    hi = m - lo + 1;
  } else {
    /* Dodge and Jureckova (2000), pag. 21. doi: 10.1007/978-1-4419-8766-2 */
    k = (m + 1) / 2. - sqrt(m);
    lo = (int) k + .5;
    k = (m + 1) / 2. + sqrt(m);
    hi = (int) k + .5;
  }
  ans *= sample_quantile(resid, m, hi) - sample_quantile(resid, m, lo);

  *value = ans;
}
