#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* These routines control the real-space one-halo term.
 * For specifics, see:
 *
 * Berlind, A.\ A., \& Weinberg, D.\ H.\ 2002, \apj, 575, 587
 * Zheng, Z. 2003, \apj, 610, 61
 * Tinker, Weinberg, Zheng, Zehavi astro-ph/0411777 (App B)
 *
 */

/* Local functions.
 */
void calc_galaxy_matter_one_halo(double *r, double *xi, int n);
double func1_xgm(double m);
double func_proj_xigm(double z);

double *xi_cs_gm2,*xi_ss_gm2,*xi_rad_gm2;

/* These are the local globals to use during the qromo integration
 */
double r_gm2, r_ds1;


/* This function tabulates the one-halo real-space term for spline interpolation.
 * If the requested radius is out-of-bounds of the tabulated function, a value of
 * zero is returned.
 */
double one_halo_galaxy_matter(double r)
{
  static int flag=0;
  static double *x,*y,*y2;
  int i,n=100;
  double a;

  if(!HOD.pdfs)return(0);

  if(!flag || RESET_FLAG_XGM1)
    {
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      RESET_FLAG_XGM1=0;
      calc_galaxy_matter_one_halo(x,y,n);
      spline(x,y,n,2.0E+30,2.0E+30,y2);
    }
  if(r>x[n])return(0);
  if(r<x[1])return(0);
  splint(x,y,y2,n,r,&a);
  return(a);

}

/* Here we calculate the one-halo real space term 
 * logarithmically spaced in r. The minimum value of r = 0.01 Mpc/h. The maximum
 * value of r is set to be approximately twice the virial radius of M_max.
 *
 * Only halos with virial radii greater than 1/2 the separation
 * contribute to the 1-halo term. 
 * Terminate integrations when r>2*R_vir(M_max).
 */
void calc_galaxy_matter_one_halo(double *r, double *xi, int n)
{
  static int ncnt=0;
  double fac,s1,rhi=1,rlo=-2,dr,mlo,x1,x2;
  int i,j;
  FILE *fp;
  char fname[100];

  ncnt++;
  rlo=log(0.001);
  rhi=log(1.9*pow(3*HOD.M_max/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0));
  dr=(rhi-rlo)/(n-1);

  if(OUTPUT>1)
    printf("calc_one_halo_gm> starting...\n");

  for(i=1;i<=n;++i)
    {
      r[i] = exp((i-1)*dr + rlo);
      xi[i] = 0;
    }
  for(i=1;i<=n;++i)
    {
      r_gm2=r[i]=exp((i-1)*dr + rlo);
      fac=1.0/(4*PI*r_gm2*r_gm2*GALAXY_DENSITY);

      mlo = 4./3.*PI*RHO_CRIT*DELTA_HALO*OMEGA_M*pow(r[i]*.5,3.0);
      if(mlo<HOD.M_low) 
	mlo = HOD.M_low;

      s1=fac*qromo(func1_xgm,log(mlo),log(HOD.M_max),midpnt);

      xi[i]=s1;
      if(OUTPUT>1)
	printf("calc_one_halo_gm> %f %e %e\n",r[i],s1,fac);
      if(s1==0) ERROR_FLAG = 0;
      if(s1<1.0E-10)break;
    }
}

/* This is the function passed to qromo in the above routine. 
 * It is the number density of
 * galaxy pairs in halos of mass m at separation r_gm2.
 * See Equation (11) from Berlind & Weinberg.
 */
double func1_xgm(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat,ss;

  if(isnan(m))
    printf("NAN func1_xgm %e\n",r_gm2);

  m=exp(m);
  cvir=halo_concentration(m)*CVIR_FAC;

  n=dndM_interp(m);
  
  nsat=N_sat(m);
  ncen=N_cen(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),THIRD);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  f_ss=dFdx_ss(r_gm2/rvir,cvir)*nsat;
  f_cs=dFdx_cs(r_gm2/rvir,cvir)*ncen;
  x=n*(f_ss+f_cs)/rvir*m/(RHO_CRIT*OMEGA_M)*m;
  return(x);

}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/* Project the galaxy-matter cross-correlation function
 * in order to acheive Delta Sigma.
 */
double delta_sigma(double rr)
{
  static int flag=0;
  double fac,s1,dr,mlo,x1,x2,mass;
  int i,j;  
  double rlo, rhi, dlogr, rmax = 20, dlogr2, a, rlo_all;

  static double *x,*y,*y2,*r,*z,*z2;
  static int n = 100;

  // first project xi_gm along the line of sight

  if(!flag || RESET_FLAG_DS)
    {
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	  r=dvector(1,n);
	  z=dvector(1,n);
	  z2=dvector(1,n);
	}
      flag=1;
      RESET_FLAG_DS=0;

      // give a call to xigm to tabulate them
      RESET_FLAG_XGM1 = 1;
      one_halo_galaxy_matter(0.1);
      if(wpx.calculate_two_halo)
	{
	  RESET_FLAG_2H = 1;
	  EXCLUSION = 5;
	  two_halo_real_space(0.1);
	  EXCLUSION = 4;
	}

      //set up limits for integration
      rlo_all = rlo = 0.002;
      rhi = 10.0;
      dlogr = log(rhi/rlo)/(n-1);

      mass = pow(10.0,0.5*(wpl.mstar_upper+wpl.mstar_lower));

      // tabulate projected xi_gm
      for(i=1;i<=n;++i)
	{
	  x[i] = exp(dlogr*(i-1))*rlo;
	  r_ds1 = x[i];
	  y[i] = 2*qromo(func_proj_xigm,0.0,rmax,midpnt);
	  //printf("tab_xigm_progj> %e %e\n",x[i],y[i]);
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);

      // integrate the projected xi_gm to get mean <R
      rlo = 0.008;
      rhi = 50.0;
      dlogr = log(rhi/rlo)/(n-1);
      for(i=1;i<=n;++i)
	{
	  r[i] = exp(dlogr*(i-1))*rlo;
	  dlogr2 = log(r[i]/rlo_all)/n;
	  s1 = 0;
	  for(j=1;j<=n;++j)
	    {
	      x1 = exp((j-0.5)*dlogr2)*rlo_all;
	      splint(x,y,y2,n,x1,&a);
	      s1 += x1*a*x1*dlogr2;
	      //printf("tab_xigm_avg> %d %e %e %e\n",i,x1,a,s1);
	    }
	  splint(x,y,y2,n,r[i],&a);
	  z[i] = 2/(r[i]*r[i]-rlo_all*rlo_all)*s1 - a;
	  z[i] = z[i]*RHO_CRIT*OMEGA_M;

	  // add in the point source
	  z[i] += mass/(r[i]*r[i]*PI);

	  // convert to solar masses /physical-pc^2
	  // r_ph = r_co*a = r_co/(1+z)
	  z[i] = z[i]*HUBBLE*1.0E-12*(1+REDSHIFT)*(1+REDSHIFT);
	  
	  //printf("tab_xigm_int> %e %e %e %e\n",r[i],2/(r[i]*r[i]-rlo_all*rlo_all)*s1,
	  //	 a,one_halo_galaxy_matter(r[i]));
	  //exit(0);
	}
      spline(r,z,n,2.0E+30,2.0E+30,z2);	 
	  
    }
  if(rr>r[n])return(0);
  if(rr<r[1])return(0);
  splint(r,z,z2,n,rr,&a);
  return(a);

}

double func_proj_xigm(double z)
{
  double rr, x;
  rr = sqrt(r_ds1*r_ds1 + z*z);
  x = one_halo_galaxy_matter(rr);
  if(wpx.calculate_two_halo)
    x += two_halo_real_space(rr);
  return x;
}

