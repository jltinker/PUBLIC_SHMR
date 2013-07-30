#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

int deconvolve()
{
  int i,j,k,n,m,npad;
  double dm, mlo, mhi, sigma = 0.15;
  double *mass, *smf, *rfunc, *decon, a[10], xx;

  a[1]=  7.455350e-02;
  a[2]=  1.188741e+01;
  a[3]=  1.739396e+00;
  a[4]=  6.541503e-01;
  a[5]=  -5.435973e+00;

  for(i=1;i<=n;++i)
    {
      rfunc[i] = 0;
      mass[i] = mlo + (i-1)*dm;
      // use the analytic model fit
      xx = pow(10.0,(mass[i]-a[2]));
      smf[i] = a[1]*pow(xx,a[3])*pow(1+pow(xx,a[4]),(a[5]-a[3])/a[4])*log(10);

      if(i>n-npad)smf[i] = 0;
    }


}
