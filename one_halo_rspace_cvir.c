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
void calc_real_space_one_halo(double *r, double *xi, int n);
double func1(double m);
double func1a(double m);
double func1b(double m);
double func1_xcorr(double m);


/* These are the local globals to use during the qromo integration
 */
double r_g2,cvir_g2,rvir_g2,slogc_g2;


/* This function tabulates the one-halo real-space term for spline interpolation.
 * If the requested radius is out-of-bounds of the tabulated function, a value of
 * zero is returned.
 */
double one_halo_real_space(double r)
{
  static int flag=0;
  static double *x,*y,*y2;
  int i,n=100;
  double a;

  if(!HOD.pdfs)return(0);

  if(!flag || RESET_FLAG_1H)
    {
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      RESET_FLAG_1H=0;
      calc_real_space_one_halo(x,y,n);
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
void calc_real_space_one_halo(double *r, double *xi, int n)
{
  double fac,s1,rhi=1,rlo=-2,dr,mlo;
  int i,j;

  rlo=log(0.01);
  rhi=log(1.9*pow(3*HOD.M_max/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0));
  dr=(rhi-rlo)/(n-1);
  
  if(OUTPUT>1)
    printf("calc_one_halo> starting...\n");
  if(!XCORR)
    GALAXY_DENSITY2 = GALAXY_DENSITY;

  for(i=1;i<=n;++i)
    {
      r_g2=r[i]=exp((i-1)*dr + rlo);
      fac=1.0/(2*PI*r_g2*r_g2*GALAXY_DENSITY*GALAXY_DENSITY2);

      mlo = 4./3.*PI*RHO_CRIT*DELTA_HALO*OMEGA_M*pow(r[i]*.5,3.0);
      if(mlo<HOD.M_low)
	mlo = HOD.M_low;

      if(XCORR)
	s1=fac*qromo(func1_xcorr,log(mlo),log(HOD.M_max),midpnt)*0.5;
      else
	s1=fac*qromo(func1,log(mlo),log(HOD.M_max),midpnt);
	
      xi[i]=s1;
      if(OUTPUT>1)
	printf("calc_one_halo> %f %e %e\n",r[i],s1,fac);
    }
}

/* This is the function passed to qromo in the above routine. 
 * It is the number density of
 * galaxy pairs in halos of mass m at separation r_g2.
 * See Equation (11) from Berlind & Weinberg.
 */
double func1(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat,slogc;

  m=exp(m);
  cvir_g2 = cvir=halo_concentration(m)*CVIR_FAC;
  slogc_g2 = slogc = 0.13*log(10);

  n=dndM_interp(m);
  
  nsat=N_sat(m);
  ncen=N_cen(m);
  

  rvir_g2 = rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  f_cs = qromo(func1a,log(cvir)-3*slogc,log(cvir)+3*slogc,midpnt)*moment_ss(m)*0.5;
  f_ss = qromo(func1b,log(cvir)-3*slogc,log(cvir)+3*slogc,midpnt)*nsat*ncen;

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  x=n*(f_cs+f_ss)/rvir*m;

  //  if(OUTPUT==3)
  //    printf("%e %e %e %e %e\n",m,n,f_ss,f_cs,rvir);

  return(x);

}

double func1a(double c)
{
  double f_ss,pc;

  f_ss=dFdx_ss(r_g2/rvir_g2,exp(c));
  pc = 1/(sqrt(TWOPI)*slogc_g2)*exp(-(c-log(cvir_g2))*(c-log(cvir_g2))/(2*slogc_g2*slogc_g2));
  return(pc*f_ss);
}

double func1b(double c)
{
  double f_cs,pc;

  f_cs=dFdx_cs(r_g2/rvir_g2,exp(c));
  pc = 1/(sqrt(TWOPI)*slogc_g2)*exp(-(c-log(cvir_g2))*(c-log(cvir_g2))/(2*slogc_g2*slogc_g2));
  return(pc*f_cs);
}



/* Same as above, only now we're calculating the number of pairs
 * in the cross-correlation.
 *
 * NB! -- We're assuming that the satellite galaxy profiles
 * of both HODs are the same.
 */
double func1_xcorr(double m)
{
  double N,n,fac2,rvir,f_ss=0,f_cs=0,cvir,x,rfof,ncen1,nsat1,ncen2,nsat2;

  m=exp(m);
  cvir=halo_concentration(m)*CVIR_FAC;
  n=dndM_interp(m);
  
  nsat1=N_sat(m);
  ncen1=N_cen(m);

  nsat2=N_sat2(m);
  ncen2=N_cen2(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   * But for x-corr, we have c1-s2, c2-s1, s1-s2.
   */
  f_ss=dFdx_ss(r_g2/rvir,cvir)*nsat1*nsat2;
  f_cs=dFdx_cs(r_g2/rvir,cvir)*(nsat1*ncen2 + nsat2*ncen1);
  x=n*(f_ss+f_cs)/rvir*m;

  return(x);

}

