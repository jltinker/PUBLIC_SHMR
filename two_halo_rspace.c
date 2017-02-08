#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"


/* This is the two-halo term of the correlation function in both real space.
 * This is solved in Fourier space and requires the 
 * non-linear power spectrum P_nl(k), which is from the model of Smith et al.
 * (see Smith et al, astro-ph/0207664).
 * 
 * For specifics see Berlind & Weinberg (2002), Zheng (2003), Tinker et al (2005)
 *
 * Fri Apr 22 09:17:30 EDT 2005 --> Now that we have a model for the scale-dependence
 *  of the halo bias with respect to the linear matter correlation function, we
 *  have an option for using the linear P_m(k) [LINEAR_PSP==1]. See details in 
 *  halo_bias.c for the linear scale dependence.
 */

/* Internal functions.
 */
double func2(double m);
double func2_matter(double m);
double func2_xigm(double m);
double func2_sig(double sig);
double func2_cen(double m);
double func2_sat(double m);
double func5(double xk);
double func_mlimit(double m);
void calc_real_space_two_halo(double *r, double *xi, int *n);
double HOD2_two_halo_real_space(double r);

/* the restricted number density.
 */
double NG_MINUS2;

/* Globals needed for the qromo functions.
 */
double r_g1,
  k1;

/* Global checkflag to stop doing the restricted
 * number density.
 */
int checkflag;

/* This tabulates the two-halo real-space term for future 
 * spline interpolation.
 */
double two_halo_real_space(double r)
{
  static int flag=0,n=45;
  static double *x,*y,*y2;
  int i;
  double max=16,min=9,a,rvir_min;
  float x1,x2;
  FILE *fp;

  //  if(r<25)
  //  return(nbody_two_halo(r));

  if(!flag || RESET_FLAG_2H)
    {
      if(!LINEAR_PSP)
	nonlinear_sigmac(8.0);
      n=30;
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      RESET_FLAG_2H=0;
      flag=1;
      if(OUTPUT>1)
	fprintf(stdout,"Calculating real-space two-halo term...\n");
      calc_real_space_two_halo(x,y,&n);
      if(!HOD.color==2)
	check_for_smoothness(x,y,n,1.5);
      XI_MAX_RADIUS=x[n];
      if(XCORR)
	for(i=1;i<=n;++i)
	  {
	    a = HOD2_two_halo_real_space(x[i]);
	    y[i] = a*y[i];
	    if(y[i]<0)y[i]=-1;
	    else y[i] = sqrt(y[i]);
	  }
      spline(x,y,n,2.0E+30,2.0E+30,y2);
    }

  /* TEST TEST TEST*/
  //R_MIN_2HALO=0.1;
  if(r<R_MIN_2HALO)return(-1);
  if(r>XI_MAX_RADIUS)return(-1);
  //if(r<R_MIN_2HALO)return(-1);
  splint(x,y,y2,n,r,&a);
  //if(a<1.58 && r<1)return 1.58;
  

  /* This check is against the spline interpolation, which
   * sometimes has (smoothness) problems with the sharp truncation of xi_2h 
   * at small separations.
   */
  if(a<-1)return(-1);
  return(a);

}

void calc_real_space_two_halo(double *r, double *xi, int *nn)
{
  double xtemp[200],rtemp[200];
  double mlimit,s2,rlo,rhi=R_MAX_2HALO,dr,klo,tolerance=1.0e-7,s1;
  int i,j,imin=0,n;
  double t0,t1,t2,t1s=0,t2s=0;

  n=*nn;

  /* Set the minimum separation of two-halo pairs.
   */
  if(HOD.color)
    {
      while(N_avg(HOD.M_low)<0.001)
	HOD.M_low*=1.01;
    }
  
  rlo = 2.2*pow(3*HOD.M_low/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  if(EXCLUSION==3)
    rlo = rlo/1.3;
  if(EXCLUSION==4)
    rlo = rlo/2.1;
  if(EXCLUSION==5)
    rlo = rlo/2.1;

  if(rlo<0.1)rlo = 0.1;
  R_MIN_2HALO = rlo;
  
  rlo = log(rlo);
  dr=(log(rhi)-rlo)/(n-1);

  checkflag=0;
  GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

  for(i=1;i<=n;++i)
    {

      r_g1=r[i]=exp(rlo+dr*(i-1));

      if(r_g1>20)
	{
	  if(EXCLUSION==5) {
	    xi[i]=xi_interp(r_g1)*GALAXY_BIAS*bias_interp(1.0E12,r_g1)/bias_interp(1.0E12,-1);
	    continue; }
	  xi[i]=xi_interp(r_g1)*GALAXY_BIAS*GALAXY_BIAS*
	    bias_interp(1.0E12,r_g1)/bias_interp(1.0E12,-1)*
	    bias_interp(1.0E12,r_g1)/bias_interp(1.0E12,-1);
	  continue;
	}

      /* Fourier transform. If BOX_SIZE specified, then truncate
       * transform at that k. If no box, take integral to k=0.
       */
      klo = 0;
      if(BOX_SIZE)klo=1/BOX_SIZE;

      j=16;
      s1 = qromo(func5,klo,j*TWOPI/r_g1,midpnt);
      s2 = s1;
      klo = j*TWOPI/r_g1;

      while(fabs(s1)>tolerance*s2) {
	j+=16;
	s1 = qromo(func5,klo,j*TWOPI/r_g1,midpnt);
	s2 += s1;
	klo = j*TWOPI/r_g1;
      }


      /* Divide the integral by the restricted number density.
       * (This is calculated in the pk_gg_2h function.)
       */
      s2=s2/NG_MINUS2;

      /* Correction factor b/c we haven't been
       * including all halos.
       * NB-- if EXCL==5, then that's xi_gm so there's only one
       * power of the numer density.
       */
      if(EXCLUSION!=5)
	xi[i]=(NG_MINUS2/GALAXY_DENSITY/GALAXY_DENSITY)*(1+s2)-1;
      else
	xi[i]=(NG_MINUS2/GALAXY_DENSITY)*(1+s2)-1;


      if(isnan(xi[i]))xi[i]=-1;

      if(xi[i]==-1)imin=i;

      if(OUTPUT>1)
	printf("calc_2halo> %f %f %d\n",r[i],xi[i],checkflag); 
      //printf("calc_2halo> %f %f %d\n",r[i],xi[i],checkflag); 
    }  

  /* Eliminate all entries which have xi=-1
   * (to eliminate a "wavy" spline fit at small r
   */
  for(j=0,i=imin+1;i<=n;++i)
    {
      j++;
      xtemp[j]=xi[i];
      rtemp[j]=r[i];
    }
  n=j;
  for(i=1;i<=n;++i)
    {
      r[i]=rtemp[i];
      xi[i]=xtemp[i];
    }
  *nn=n;
  R_MIN_2HALO=r[1];    

  /*printf("calc_2halo> %d %d %f\n",imin,n,r[1]);*/
}


/* This calculates and tabulates the galaxy power spectrum. This is done
 * by taking the galaxy-number-weighted halo bias (with scale-dependent 
 * halo bias and extended structure of halos taken into account) and multiplying
 * it by the non-linear matter power spectrum.
 * 
 * The galaxy-number-weighted average is an integral over the halo mass function
 * (see Equations B10 & B12 in Appendix B of Tinker et al.) over all allowed halo
 * pairs.
 *
 * The allowed number of halo pairs is controlled by the type of halo exclusion used:
 *  EXCLUSION = 1 --> only halos with Rvir < r/2
 *  EXCLUSION = 2 --> only halo PAIRS with R1+R2 < r ("spherical eclusion")
 *  EXCLUSION = 3 --> a fraction of halo pairs f(x), x=r/(R1+R2). ("ellipsoidal exclusion")
 *  EXCLUSION = 4 --> only halos with Rvir < r
 *
 * 1 is fast b/c the integral is separable, calculated once and squared. Can drastically
 * underestimate the number of small-sep pairs.
 *
 * 2 and 3 are slower since it is a true double integral over halo masses M1 and M2, but
 * more accurate. Both 2 and 3 use the approximation of matching the restricted number 
 * density and using a seperable double integral (once again, see App B of Tinker et al).
 * The actual double integral for doing spherical or ellipsoidal exclusion is not at
 * present included in the public code.
 */

double psp_gg_2h(double k, double r)
{
  static double rp=-1,*x,*y,*y2;
  static int flag=0;
  int n=60,i;
  double a,dk,klo=-3,khi=3,xk,s1,s2,mhi,mlo,mhi1,mlo1,am,siglo,sighi;

  double t1,t0,ttot1,ttot2;


  /* The fourier transform is done at each r, so tabulate the power spectrum
   * at all k for a given r. If r has changed, re-tabulate.
   */
  if(rp!=r)
    {
      mlo=HOD.M_low;
      if(!flag)
	{
	  flag=1;
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      dk=(khi-klo)/n;

      switch(EXCLUSION) {
      case 1:
	mhi=4./3.*DELTA_HALO*RHO_CRIT*PI*r*r*r*OMEGA_M*0.125; 
	if(mhi>HOD.M_max)mhi=HOD.M_max;
	NG_MINUS2 = qromo(func_galaxy_density,log(mlo),log(mhi),midpnt);
	NG_MINUS2*=NG_MINUS2;
	for(i=n;i>=1;--i)
	  {
	    xk=pow(10.0,klo+dk*i);
	    k1=xk;
	    if(nfw_transform(xk,mhi)<0.999)
	      a=qromo(func2,log(mlo),log(mhi),midpnt);
	    x[i]=log(xk);
	    if(LINEAR_PSP)
	      y[i]=log(linear_power_spectrum(xk)*a*a); 
	    else
	      y[i]=log(nonlinear_power_spectrum(xk)*a*a); 
	  }
	break;
      case 2:
      case 3:
      case 4: // where max halo mass is M(R=r) rather than M(R=r/2)
	// NB: if VARIABLE_EXCLUSION=1, then this can vary between r and r/2. 
	mhi=4./3.*DELTA_HALO*RHO_CRIT*PI*r*r*r*OMEGA_M/pow(EXCLUSION_RADIUS,3.0); 
	if(mhi>HOD.M_max)mhi=HOD.M_max;
	NG_MINUS2 = qromo(func_galaxy_density,log(mlo),log(mhi),midpnt);
	NG_MINUS2*=NG_MINUS2;
	for(i=n;i>=1;--i)
	  {
	    xk=pow(10.0,klo+dk*i);
	    k1=xk;
	    if(nfw_transform(xk,mhi)<0.999)
	      {
		a=qromo(func2,log(mlo),log(mhi),midpnt);
	      }
	    //if(i==1)printf("MUH %e %e %e %e\n",r_g1,mlo,mhi,a/sqrt(NG_MINUS2));// exit(0);
	    x[i]=log(xk);
	    if(LINEAR_PSP)
	      y[i]=log(linear_power_spectrum(xk)*a*a); 
	    else
	      y[i]=log(nonlinear_power_spectrum(xk)*a*a); 
	  }
	break;

      case 5: // this is going to return the galaxy-mass cross correlation function.
	// so the minimum halo mass for the two terms is different.
	mhi=4./3.*DELTA_HALO*RHO_CRIT*PI*r*r*r*OMEGA_M; 
	if(mhi>HOD.M_max)mhi=HOD.M_max;
	NG_MINUS2 = qromo(func_galaxy_density_xigm,log(mlo),log(mhi),midpnt); //DON'T square this. (as done above)
	siglo = sigmac_interp(mlo);
	sighi = sigmac_interp(mhi);
	for(i=n;i>=1;--i)
	  {
	    xk=pow(10.0,klo+dk*i);
	    k1=xk;
	    if(nfw_transform(xk,mhi)<0.999)
	      {
		am=qromo(func2_sig,log(1/siglo)-100,log(1/sighi),midpnt);
		a=qromo(func2_xigm,log(mlo),log(mhi),midpnt);	
	      }
	    //if(i==1)printf("MUH %e %e %e %e %e\n",r_g1,mlo,mhi,am,a/NG_MINUS2);// exit(0);
	    x[i]=log(xk);
	    if(LINEAR_PSP)
	      y[i]=log(linear_power_spectrum(xk)*a*am); 
	    else
	      y[i]=log(nonlinear_power_spectrum(xk)*a*am); 
	  }
	break;

      }
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      rp=r;
    }
  k=log(k);
  splint(x,y,y2,n,k,&a);
  return(exp(a));
}


/* The routine called by qromo for the integral in Zheng's Eq [7]
 */
double func2(double m)
{
  double n,N,b,yg,x,Ncen,Nsat;

  m=exp(m);
  n=dndM_interp(m);
  b=bias_interp(m,r_g1);
  Ncen=N_cen(m);
  Nsat=N_sat(m);
  yg=nfw_transform(k1,m);
  //x=n*N*b*yg*m;
  x = n*(Ncen + Nsat*yg)*b*m;
  //printf("%e %e %f %f %e %f\n",m,n,b,x,yg,bias_interp(m,-1.0)); 
  return(x);
}

/* This routine is over all halo masses, but variable is sigma
 * rather than mass.
 */
double func2_sig(double sig)
{
  double m,n,b,yg,x;
  x = sig;
  sig = 1/exp(sig);
  n = fsigma(sig);
  b = bias_fsigma(sig,r_g1);
  m = sigma2mass(sig);
  //yg = nfw_transform(k1,m);
  yg = 1;
  //  printf("%e %e %e %e %e %f %e %f\n",x,log(1/sig),sig,m,n,b,yg,bias_fsigma(sig,-1)); 
  //if(b>100)exit(0);
  return n*b*yg;
}

/* The routine called by qromo for the integral in Zheng's Eq [7]
 */
double func2_xigm(double m)
{
  double n,N,b,yg,x,Ncen,Nsat;

  m=exp(m);
  n=dndM_interp(m);
  b=bias_interp(m,r_g1);
  Ncen=Ncen_xigm(m);
  Nsat=Nsat_xigm(m);
  yg=nfw_transform(k1,m);
  //x=n*N*b*yg*m;
  x = n*(Ncen + Nsat*yg)*b*m;
  //  printf("%e %e %f %f %e %f\n",m,n,b,N,yg,bias_interp(m,-1.0)); 
  return(x);
}

/* The routine called by qromo for the 2halo term contributed by
 * the MATTER.
 */
double func2_matter(double m)
{
  double n,N,b,yg,x,Ncen,Nsat;

  m=exp(m);
  //n=dndM_interp(m);
  fmuh(m);
  n = halo_mass_function(m);
  b=bias_interp(m,r_g1);
  yg=nfw_transform(k1,m);
  x=n*b*yg*m*m/(OMEGA_M*RHO_CRIT);
  //printf("%e %e %f %e %f\n",m,n,b,yg,bias_interp(m,-1.0)); 
  return(x);
}


/* The routine called by qromo for the integral in Zheng's Eq [7]
 * FOR CENTRAL GALAXIES ONLY: This means no Fourier transform of NFW profile.
 */
double func2_cen(double m)
{
  double n,N,b,yg,x;

  m=exp(m);
  n=dndM_interp(m);
  b=bias_interp(m,r_g1);
  N=N_cen(m);
  x=n*N*b*m;
  return(x);
}

/* The routine called by qromo for the integral in Zheng's Eq [7]
 * FOR SATELLITE GALAXIES ONLY.
 */
double func2_sat(double m)
{
  double n,N,b,yg,x;

  m=exp(m);
  n=dndM_interp(m);
  b=bias_interp(m,r_g1);
  N=N_sat(m);
  yg=nfw_transform(k1,m);
  x=n*N*b*yg*m;
  return(x);
}


/* This is the integrand of the Fourier transform
 * of the two-halo power spectrum to correlation function
 */
double func5(double xk)
{
  double psp1;
  double xk1,xk2,x3;

  if(xk==0)return(0);
  psp1=psp_gg_2h(xk,r_g1);
  xk1=r_g1*xk;
  psp1*=sin(xk1)/xk1/xk;
  return(psp1);
}

/* This is the function sent to zbrent to calculate the halo mass 
 * which gives the matched restricted number density.
 */
double func_mlimit(double m)
{
  double s1;
  if(N_avg(exp(m))>0)
    s1=qromo(func_galaxy_density,log(HOD.M_low),m,midpnt);
  else
    s1=0;
  return(s1-NG_MINUS2);
}


/* This tabulates the two-halo real-space term
 * for the second HOD function for future 
 * spline interpolation.
 */
double HOD2_two_halo_real_space(double r)
{
  static int flag=0,n=45;
  static double *x,*y,*y2;
  int i;
  double max=16,min=9,a,rvir_min,galtemp;
  float x1,x2;
  FILE *fp;

  if(!flag || RESET_FLAG_2H)
    {
      /* Switch HODs temporarily
       */
      HODt.M_min = HOD.M_min;
      HODt.M_low = HOD.M_low;
      HODt.M1 = HOD.M1;
      HODt.alpha = HOD.alpha;
      HODt.M_cen_max = HOD.M_cen_max;
      HODt.sigma_logM = HOD.sigma_logM;
      HODt.M_max = HOD.M_max;
      HODt.pdfc = HOD.pdfc;
      HODt.pdfs = HOD.pdfs;
      galtemp = GALAXY_DENSITY;
      
      HOD.M_min = HOD2.M_min;
      HOD.M_low = HOD2.M_low;
      HOD.M1 = HOD2.M1;
      HOD.alpha = HOD2.alpha;
      HOD.M_cen_max = HOD2.M_cen_max;
      HOD.sigma_logM = HOD2.sigma_logM;
      HOD.M_max = HOD2.M_max;
      HOD.pdfc = HOD2.pdfc;
      HOD.pdfs = HOD2.pdfs;
      GALAXY_DENSITY = GALAXY_DENSITY2;
      set_HOD_params();

      if(!LINEAR_PSP)
	nonlinear_sigmac(8.0);
      n=30;
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      RESET_FLAG_2H=0;
      flag=1;
      if(OUTPUT)
	fprintf(stdout,"Calculating HOD2 real-space two-halo term...\n");
      calc_real_space_two_halo(x,y,&n);
      check_for_smoothness(x,y,n,1.5);
      XI_MAX_RADIUS=x[n];
      spline(x,y,n,2.0E+30,2.0E+30,y2);

      /* Switch HODs back
       */
      HOD.M_min = HODt.M_min;
      HOD.M_low = HODt.M_low;
      HOD.M1 = HODt.M1;
      HOD.alpha = HODt.alpha;
      HOD.M_cen_max = HODt.M_cen_max;
      HOD.sigma_logM = HODt.sigma_logM;
      HOD.M_max = HODt.M_max;
      HOD.pdfc = HODt.pdfc;
      HOD.pdfs = HODt.pdfs;
      GALAXY_DENSITY = galtemp;
      set_HOD_params();
    }

  if(r>XI_MAX_RADIUS)return(-1);
  if(r<R_MIN_2HALO)return(-1);
  splint(x,y,y2,n,r,&a);

  /* This check is against the spline interpolation, which
   * sometimes has (smoothness) problems with the sharp truncation of xi_2h 
   * at small separations.
   */
  if(a<-1)return(-1);
  return(a);

}

