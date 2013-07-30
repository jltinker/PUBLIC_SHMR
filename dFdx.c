#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "header.h"

/** halo pair separation profile for halo with concentration c_NFW
    [satellite-satellite]
    x=r/(2*R_vir)  0<=x<=1
    analytic solution of Sheth et al (2001) MNRAS 325, 1288 (eq.A25)
    dF/dx \propto \lambda(r)*r^2
    normalized so that \int_0^1 dF/dx dx =1
    The normalization factor in this subroutine is obtained from my 
    fitting formula, which has a fractional error less than 0.2% for
    c_NFW in the range ~1 to ~100. The formula is basically a double
    powerlaw with a smooth transition and the residuals are further
    reduced with a (1+sine) function.
 **/
 
double dFdx_ss(double x, double c_NFW)
{
 double f,y,a,Anorm;
 double A0=4.915,alpha=-3.099,beta=0.617,cNFWc=1.651,mu=4.706;
 double B0=0.0336,omega=2.684,phi=0.4079;
 /** above are parameters in fitting formula of the normalization **/
 double t;
 double t1,t2,t3,ta;

 if(x<=0.0 || x>=1.0) return 0.0;

 a=1.0/c_NFW;
 y=2.0*c_NFW*x;


 if(x<0.5)
   {
     f=(-4.0*(1.0+a)+2.0*a*y*(1.0+2.0*a)+a*a*y*y)/(2.0*(1.0+a)*(1.0+a)*(2.0+y));
     f=f+log(fabs((1.0+a-a*y)*(1.0+y))/(1.0+a))/y+y*log(1.0+y)/(2.0+y)/(2.0+y);
   }
 else
   {
     f=y*log((1.0+a)/fabs(a*y+a-1))/(2.0+y)/(2.0+y);
     f=f+0.5*y*(a*a*y-2.0*a)/(1.0+a)/(1.0+a)/(2.0+y);
   }

 /** get the normalization factor from a fitting formula **/
 t=pow(c_NFW/cNFWc,(beta-alpha)/mu);
 Anorm=A0*pow(c_NFW,alpha)*pow(1.0+t,mu)*(1.0+B0*sin(omega*(log10(c_NFW)-phi)));

 return Anorm*f;
}

/** [central-satellite] i.e. the NFW profile rho(r)*r^2
    x=r/(2*R_vir)
 **/

double dFdx_cs(double x, double c_NFW)
{
 double f,A,y;

 if(x>0.5) return 0.0;
 else {
   y=1.0+2.0*c_NFW*x;  // (1+r/rs)
   f=x/(y*y);
   A=(log(1.0+c_NFW)-c_NFW/(1.0+c_NFW))/(4.0*c_NFW*c_NFW);
   return f/A;
 }
}

/* The redshift-space one-halo term spends the majority of its time 
 * calculating NFW pair densities for sat-sat pairs. To speed up
 * the calculation, this tabulates the values of dFdx_ss on a grid
 * for bilinear interpolation.
 *
 * In tests, accuracy is usually about 1E-5, and at worst 1E-3.
 */
double dFdx_ss_interp(double r, double c)
{
  static int flag=0,reset=0;
  static double **x;
  static double clo,dlogc,chi;
  int nx=2000, nc=100, i,j,ir,ic;
  double c_nfw,x1,x2,x3,c_fac=2;

  if(!flag)
    {
      x=dmatrix(0,nc,1,nx);
      flag=1;
      if(OUTPUT)
	fprintf(stdout,"dFdx_ss_interp> Tabulating dFdx_ss...\n");
      clo = 0.1;
      chi = 1000.0;
      dlogc = (log(chi)-log(clo))/(nc-1);
      for(j=1;j<=nx;++j)
	x[0][j]=0;
      for(i=1;i<=nc;++i)
	{
	  c_nfw=exp((i-1)*dlogc)*clo;
	  for(j=1;j<=nx;++j) { 
	    x[i][j]=dFdx_ss((double)j/nx,c_nfw); }
	}
      if(OUTPUT)
	fprintf(stdout,"dFdx_ss_interp> ...Finished.\n");
    }
  r*=nx;
  c_fac = log(c/clo)/dlogc+1;
  ir=(int)r;
  ic=(int)c_fac;
  if(ic==0)ic++;
  if(ic==nc)ic--;
  if(ic==0 || ic>=nc) {
    printf("%f %f %f %d\n",r/nx,c,c_fac,ic);
    endrun("dFdx error"); }

  x1=(x[ic][ir+1]-x[ic][ir])*(r-ir)+x[ic][ir];
  x2=(x[ic+1][ir+1]-x[ic+1][ir])*(r-ir)+x[ic+1][ir];
  return((x2-x1)*(c_fac-ic)+x1);
}


/* Well, it looks like it really doesn't help much of anything to tabulate
 * the central-satellite function (which isn't too surprising since the number
 * of computations isn't very large). But the above one seems to do very well.
 */

double dFdx_cs_interp(double r, double c)
{
  static int flag=0;
  static double **x;
  int nx=1000, nc=19, i,j,ir,ic;
  double c_nfw,x1,x2,x3;

  if(r>0.5)return(0);

  if(!flag++)
    {
      fprintf(stderr,"Tabulating dFdx_cs...\n");
      x=dmatrix(0,nc,1,nx);
      for(i=1;i<=nc;++i)
	{
	  c_nfw=i;
	  x[0][j]=0;
	  for(j=1;j<=nx;++j)
	    x[i][j]=dFdx_cs((double)j/nx,c_nfw);
	  fprintf(stderr,"%d\n",i);
	}
      fprintf(stderr,"Finished.\n");
    }
  r*=nx;
  ir=(int)r;
  ic=(int)c;
  if(ic==0 || ic>=nc)
    endrun("dFdx error");

  x1=(x[ic][ir+1]-x[ic][ir])*(r-ir)+x[ic][ir];
  x2=(x[ic+1][ir+1]-x[ic+1][ir])*(r-ir)+x[ic+1][ir];
  return((x2-x1)*(c-ic)+x1);
}
