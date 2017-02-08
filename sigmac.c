#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/*
 *    sigmac calls qromo to evaluate the integral

         int_0^infty P_s(k) 4 pi k^2 dk / (2pi)^3

      where P_s(k) is the power spectrum smoothed through a gaussian
      window of radius rgaus and a top hat window of radius rth.
      The power spectrum is a power law of index xindx, multiplied by
      the square of the CDM transfer function if itrans = 1.

      If the smoothing radius [rad] is positive, a top-hat window
      function is used. 
      If [rad<0], a gaussian window function is used.

*/


double rad1;
double *g22_x, *g22_y, *g22_y2, g22_sig;
int g22_n;
double func_sigmac(double xk);
double func_sigmac_nl(double xk);
double sigmac_sirko(double rad);


double func_s2m(double m)
{
  double a;
  splint(g22_x,g22_y,g22_y2,g22_n,m,&a);
  return a-g22_sig;
}
double sigma2mass(double sigma)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2, pnorm;
  int i,n=100;
  static double dlogm,max=5.e17,min=1.0e-7,a,b,m1,m2,dm1,xi,power,rm,sig,b1,b2,mass;
  double m;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  g22_x=dvector(1,n);
	  g22_y=dvector(1,n);
	  g22_y2=dvector(1,n);
	}
      flag=1;
      dlogm = (log(max) - log(min))/(n-1);
      pnorm=SIGMA_8/sigmac(8.0);
      for(i=1;i<=n;++i)
	{
	  m = exp((i-1)*dlogm)*min;
	  rm=pow(3.0*m/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
	  sig=log(1/(pnorm*sigmac(rm)));
	  g22_x[i] = log(m);
	  g22_y[i] = sig;
	  //printf("SIG %e %e\n",exp(g22_x[i]),exp(g22_y[i]));
	}
      g22_n = n;
      spline(g22_x,g22_y,g22_n,2.0E+30,2.0E+30,g22_y2);
      prev_cosmo=RESET_COSMOLOGY;

      min = log(min);
      max = log(max);
    }
  if(log(1/sigma)<g22_y[1])return exp(g22_x[1]);
  if(log(1/sigma)>g22_y[n]) {
    printf("BOO %e %e\n",sigma,exp(g22_y[n])); }
  g22_sig = log(1/sigma);
  m = zbrent(func_s2m,min,max,1.0E-5);
  return exp(m);
}

double sigmac_interp(double m)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2, pnorm;
  int i,n=100;
  double dlogm,max=5.e16,min=1.0e7,a,b,m1,m2,dm1,xi,power,rm,sig,b1,b2,mass;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dlogm = (log(max) - log(min))/(n-1);
      pnorm=SIGMA_8/sigmac(8.0);
      for(i=1;i<=n;++i)
	{
	  mass = exp((i-1)*dlogm)*min;
	  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
	  sig=log(pnorm*sigmac(rm));
	  x[i] = log(mass);
	  y[i] = sig;
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;

    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return exp(a);
}

double sigmac_radius_interp(double m)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2, pnorm;
  int i,n=100;
  double dlogm,max=80.0,min=0.1,a,b,m1,m2,dm1,xi,power,rm,sig,b1,b2,mass;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dlogm = (log(max) - log(min))/(n-1);
      pnorm=SIGMA_8/sigmac(8.0);
      for(i=1;i<=n;++i)
	{
	  rm = exp((i-1)*dlogm)*min;
	  sig=log(pnorm*sigmac(rm));
	  x[i] = log(rm);
	  y[i] = sig;
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;

    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return exp(a);
}


double sigmac(double rad)
{
  double xk1,s1=1,s2=0,sigma0,lowL=0;

  if(HIGH_PRECISION)
    return sigmac_sirko(rad);

  rad1=rad;
  if(BOX_SIZE>0)lowL = TWOPI/BOX_SIZE;
  xk1 = 1./(2.*PI*mabs(rad)); 
  s1=qromo(func_sigmac,lowL,xk1,midpnt);
  s2=qromo(func_sigmac,xk1,1.e20,midinf);
  sigma0=sqrt((s1+s2)*(4*PI)/pow(2*PI,3.0));
  return sigma0; 
  
}

double sigmac_sirko(double rad)
{
  double xk1,s1=1,s2=0,sigma0;
  double integ_range,total,answer,allowed_error=1.0e-10;
  int izone;

  integ_range = 1.*PI/rad;
  answer = 1.; // dummy value to force while loop
  total = 0.;
  izone = 0;
  while (fabs(answer) > allowed_error*fabs(total)) {
    answer = integrate(func_sigmac,izone*integ_range,(izone+1)*integ_range,1,1);
    //answer = qromo(func_sigmac,izone*integ_range,(izone+1)*integ_range,midpnt);
    total += answer;
    izone++;
  }
  return(sqrt(total*(4*PI)/pow(2*PI,3.0)));
}



double func_sigmac(double xk)
{
  double xkr,w,psp;

  if(rad1>0)
    {
      xkr = xk*rad1;
      w = 3*(sin(xkr)-xkr*cos(xkr))/(xkr*xkr*xkr);
      w = w*w;
    }
  else
    {
      xkr = -rad1*xk;
      w = exp(-xkr*xkr);
    }

  if(ITRANS>0)
    psp=pow(xk,SPECTRAL_INDX + SPECTRAL_RUN*log(xk/0.002))*pow(transfnc(xk),2.0);
  else
    psp=pow(xk,SPECTRAL_INDX + SPECTRAL_RUN*log(xk/0.002));
  psp=psp*w*xk*xk;
  return(psp);
}

/* Same as above, but no instead of using the linear
 * theory power spectrum, use the non-linear P(k)
 */
double nonlinear_sigmac(double rad)
{
  double xk1,s1=1,s2=0,sigma0;

  rad1=rad;
  xk1 = 1./(2.*PI*mabs(rad)); 
  s1=qromo(func_sigmac_nl,0.0,xk1,midpnt);
  s2=qromo(func_sigmac_nl,xk1,1.e20,midinf);
  sigma0=sqrt((s1+s2));
  return sigma0; 
}

double func_sigmac_nl(double xk)
{
  double xkr,w,psp;

  if(rad1>0)
    {
      xkr = xk*rad1;
      w = 3*(sin(xkr)-xkr*cos(xkr))/(xkr*xkr*xkr);
      w = w*w;
    }
  else
    {
      xkr = -rad1*xk;
      w = exp(-xkr*xkr);
    }

  psp=nonlinear_power_spectrum(xk);
  psp=psp*w/xk;
  return(psp);
}
