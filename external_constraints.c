#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

double external_constraints()
{
  double mlo, mhi, chi2,chi2s;
  static int niter=0;
  /* What is the bias of this model at logM=9.5?
   */
  mlo = pow(10.0,9.9);
  mhi = pow(10.0,10.1);
  set_up_hod_for_shmr(mlo,mhi,wpl.a);
  
  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

  chi2 = (GALAXY_BIAS-1.2)*(GALAXY_BIAS-1.2)/(0.12*0.12);

  if(EXTERNAL_CONSTRAINTS>1)
    {
      chi2 = 0;
      chi2s = (wpl.a[6]-0.245)*(wpl.a[6]-0.245)/(0.001*0.001);
      chi2 += chi2s;
      chi2s = (wpl.a[3]-0.38)*(wpl.a[3]-0.38)/(0.05*0.05);
      chi2 += chi2s;
    }
  printf("EXTERN%d %e\n",niter++,chi2);
  return chi2;
}
