#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

double collapse_redshift(double z);
double cvir_pnorm_g1,
  cvir_sigma_g1;

/* This calculates and tabulates the halo concentrations
 * as a function of halo mass. Uses the "Bullock model", 
 * described in a little more detail below.
 */

double halo_concentration(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2,cfac;
  double a,dm,x3,x4;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      MSTAR = mstar();
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",DELTA_HALO);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i] = log(12*pow(x[i]/1.0e12,-0.11)); //zentner comparison neto et al paper.
	  x[i] = log(x[i]);
	  continue;
	  //y[i]=munoz_cuartas_cvir_model(x[i]);
	  // convert from cvir to c_delta
	  x[i]=log(halo_mass_conversion(x[i],&y[i],DELTA_HALO));
	  y[i]=log(y[i]);
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));
}


/* Model from Munoz_cuartas 2010 , WMAP5, Cvir mass = h^-1
 */
double munoz_cuartas_cvir_model(double mass)
{
  double alpha, beta, gamma, b_cvir, a_cvir, w_cvir, m_cvir,cvir;
  
  w_cvir = 0.029;
  m_cvir = 0.097;
  alpha  = -110.001;
  beta   = 2469.720;
  gamma  = 16.885;

  b_cvir = (alpha/(REDSHIFT+gamma))+(beta/pow(REDSHIFT+gamma,2));
  a_cvir = w_cvir*REDSHIFT - m_cvir;

  cvir = pow(10.0, a_cvir*log10(mass)+b_cvir);
  //fprintf("cvir and mass %e %e \n",cvir, mass);

  return(cvir);
}


/* Some quantities are specific to an overdensity of 200 (i.e., the Jenkins mass
 * function and the halo bias relation of Tinker et al. )
 * Specfically, these are actually for FOF 0.2 halos, for which 200 is the
 * current best approximation. (Although Warren et al quotes 250, which is the most recent
 * value.)
 *
 * Therefore, the halo concentrations for Delta=200 halos need to be easily 
 * accessible for halo mass conversion of these quantities. The values are 
 * tabulates here, as opposed to the above routine which tabulates the 
 * concentrations for a user-specified overdensity.
 */

double halo_c200(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2;
  double a,dm,x3,x4;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",200.0);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=munoz_cuartas_cvir_model(x[i]);
	  x[i]=log(halo_mass_conversion(x[i],&y[i],200.0));
	  y[i]=log(y[i]);
	}
    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));
  
}
