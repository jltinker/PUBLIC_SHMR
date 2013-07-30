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
double func1cs(double m);
double func1satsat(double m);
double one_halo_ss(double r);
double one_halo_cs(double r);

double *xi_cs_g2,*xi_ss_g2,*xi_rad_g2;

/* These are the local globals to use during the qromo integration
 */
double r_g2;


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
  static int ncnt=0;
  double fac,s1,rhi=1,rlo=-2,dr,mlo,x1,x2;
  int i,j;
  FILE *fp;
  char fname[100];

  ncnt++;
  rlo=log(0.01);
  rhi=log(1.9*pow(3*HOD.M_max/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0));
  dr=(rhi-rlo)/(n-1);
  
  if(OUTPUT>1)
    printf("calc_one_halo> starting...\n");
  GALAXY_DENSITY2 = GALAXY_DENSITY;

  if(OUTPUT>2)
    {
      sprintf(fname,"%s.1halo",Task.root_filename);
      fp = fopen(fname,"w");
    }

  for(i=1;i<=n;++i)
    {
      r[i] = exp((i-1)*dr + rlo);
      xi[i] = 0;
    }
  for(i=1;i<=n;++i)
    {
      r_g2=r[i]=exp((i-1)*dr + rlo);
      fac=1.0/(2*PI*r_g2*r_g2*GALAXY_DENSITY*GALAXY_DENSITY2);

      mlo = 4./3.*PI*RHO_CRIT*DELTA_HALO*OMEGA_M*pow(r[i]*.5,3.0);
      if(mlo<HOD.M_low)
	mlo = HOD.M_low;

      s1=fac*qromo(func1,log(mlo),log(HOD.M_max),midpnt);

      xi[i]=s1;
      if(OUTPUT>1)
	printf("calc_one_halo> %f %e %e\n",r[i],s1,fac);
      if(s1==0) ERROR_FLAG = 0;
      if(s1<1.0E-10)break;
      continue;

      x1 = fac*qromo(func1satsat,log(mlo),log(HOD.M_max),midpnt);
      x2=0;
      if(r_g2<exp(rhi)/2)
	x2 = fac*qromo(func1cs,log(mlo),log(HOD.M_max),midpnt);
      if(OUTPUT>1)printf("MOO%d %f %e %e\n",ncnt,r_g2,x1,x2);
      if(OUTPUT>2)fprintf(fp,"%f %e %e\n",r_g2,x1,x2);

    }
  if(OUTPUT>2)fclose(fp);
}

/* This is the function passed to qromo in the above routine. 
 * It is the number density of
 * galaxy pairs in halos of mass m at separation r_g2.
 * See Equation (11) from Berlind & Weinberg.
 */
double func1(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat,ss;

  m=exp(m);
  cvir=halo_concentration(m)*CVIR_FAC;

  n=dndM_interp(m);
  
  nsat=N_sat(m);
  ncen=N_cen(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  f_ss=dFdx_ss(r_g2/rvir,cvir)*moment_ss(m)*0.5;
  f_cs=dFdx_cs(r_g2/rvir,cvir)*nsat*ncen;
  x=n*(f_ss+f_cs)/rvir*m;
  return(x);
}

double func1satsat(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat;

  m=exp(m);
  cvir=halo_concentration(m)*CVIR_FAC;

  n=dndM_interp(m);
  
  nsat=N_sat(m);
  ncen=N_cen(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  f_ss=dFdx_ss(r_g2/rvir,cvir)*moment_ss(m)*0.5;
  x=n*(f_ss)/rvir*m;

  //  if(OUTPUT==3)
  //    printf("%e %e %e %e %e\n",m,n,f_ss,f_cs,rvir);

  return(x);

}

double func1cs(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat;

  m=exp(m);
  cvir=halo_concentration(m)*CVIR_FAC;

  n=dndM_interp(m);
  
  nsat=N_sat(m);
  ncen=N_cen(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  f_cs=dFdx_cs(r_g2/rvir,cvir)*nsat*ncen;
  x=n*(f_cs)/rvir*m;

  //  if(OUTPUT==3)
  //    printf("%e %e %e %e %e\n",m,n,f_ss,f_cs,rvir);

  return(x);

}


