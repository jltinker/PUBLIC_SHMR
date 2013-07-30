#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

void mcmc_clustering_input(void);
double biasfit(int isample, double bb);
double find_unbiased_hod(double mmin);

void fit_for_bias()
{
  int i,j,k, nb=1000, ichi2min;
  double chi2, xbias[nb], xchi2[nb], chi2min;
  double blo, bhi, db, om, xm, error_lo, error_hi;

  mcmc_clustering_input();

  // find an HOD with bias=1
  HOD.M1 = 1.0E+13;
  HOD.M_cut = 1.0E+14;
  HOD.alpha = 1;
  HOD.pdfc = 1;
  HOD.pdfs = 1;
  zbrent(find_unbiased_hod,log(1.0E+10),log(1.0E+13),1.0E-3);
  fprintf(stdout,"unbiased HOD: %e\n",HOD.M_min);
  blo = 0.2;
  bhi = 10.0;
  db = (bhi-blo)/(nb-1);
  HOD.M_low = HOD.M_min;

  om = (OMEGA_M*pow(1+REDSHIFT,3))/(OMEGA_M*pow(1+REDSHIFT,3) + (1-OMEGA_M));
  BETA = pow(om,GAMMA)/GALAXY_BIAS;
  printf("BETA: %f %f\n",om,BETA);
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;
  one_halo_real_space(1) + two_halo_real_space(1);

  for(i=1;i<=wpl.nsamples;++i)
    {
      chi2min = 1.0E+6;
      for(j=0;j<nb;++j)
	{
	  xbias[j] = blo + j*db;
	  xchi2[j] = biasfit(i,xbias[j]);
	  if(xchi2[j]<chi2min) { chi2min = xchi2[j]; ichi2min = j; }
	  //printf("%d %e %e\n",j,xbias[j],xchi2[j]);
	}
      for(j=ichi2min;j>=0;--j)
	if(xchi2[j]>chi2min+4)break;
      error_lo = xbias[ichi2min] - xbias[j];
      for(j=ichi2min;j<nb;++j)
	if(xchi2[j]>chi2min+4)break;
      error_hi = xbias[j]-xbias[ichi2min];

      printf("BIAS %d %e %e %e %e\n",i,xbias[ichi2min],error_lo, error_hi,chi2min);
      for(j=1;j<=wpl.ndata[i];++j)
	{
	  xm = xbias[ichi2min]*xbias[ichi2min]*
	    projected_xi_matter(wpl.rdata[i][j])*bias_interp(1.0E12,wpl.rdata[i][j])/bias_interp(1.0E12,-1)*
	    projected_xi(wpl.rdata[i][j])/projected_xi_rspace(wpl.rdata[i][j]);
	  printf("BIASFIT%d %e %e %e %e\n",i,wpl.rdata[i][j],wpl.xdata[i][j],wpl.edata[i][j],xm);
	}
    }
  exit(0);
}

double biasfit(int isample, double bb)
{
  int i,j;
  double xm[100], chi2 = 0, r;

  // get model: account for scale-dependent bias and RSD
  for(i=1;i<=wpl.ndata[isample];++i)
    xm[i] = bb*bb*projected_xi_matter(wpl.rdata[isample][i])*bias_interp(1.0E12,wpl.rdata[isample][i])/bias_interp(1.0E12,-1)*
      projected_xi(wpl.rdata[isample][i])/projected_xi_rspace(wpl.rdata[isample][i]);;

  for(i=1;i<=wpl.ndata[isample];++i)
    for(j=1;j<=wpl.ndata[isample];++j)
      chi2 += (xm[i]-wpl.xdata[isample][i])*(xm[j]-wpl.xdata[isample][j])*wpl.covar[isample][i][j];

  return chi2;

}

double find_unbiased_hod(double mmin)
{
  HOD.M_min = exp(mmin);
  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_min),log(HOD.M_max),midpnt);
  GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_min),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  return GALAXY_BIAS - 2.0;
}
