#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This is a calculation of the differential halo mass function of
 * Jenkins et al 2001 (MNRAS 321, 372). 
 *
 */
void initialize_mass_function2(double *a1, double *a2, double *a3, double *a4, double *a5);
void initialize_mass_function(double *a1, double *a2, double *a3, double *a4);

double every_fucking_observers_mass_function(double mass)
{
  double h = HUBBLE;
  return dndM_interp(mass/h)*pow(HUBBLE,3.0);
}

/* Just returns f(sigma) for the mass function
 */
double fsigma(double sig)
{
  static double a1, a2, a3, a4, a5;
  static double pnorm, prev_delta=-1, prev_cosmo=-1;
  if(DELTA_HALO != prev_delta || prev_cosmo != RESET_COSMOLOGY)
    {
      initialize_mass_function2(&a1,&a2,&a3,&a4,&a5);
      prev_delta = DELTA_HALO;
      prev_cosmo = RESET_COSMOLOGY;
    }
  return a5*(pow(sig/a3,-a2)+pow(sig,-a1))*exp(-a4/sig/sig);
}

double halo_mass_function(double mass)
{
  double sig,logm,a,slo,shi,rm,rlo,rhi,mlo,mhi,dsdM,n,nuprime,nufnu,p,A, fac;
  static int flag=0,SO180=0,SO324=0,WARREN=0,ST=0,JENKINS=0;
  static double pnorm, prev_delta, prev_cosmo;
  double btemp = -1;

  static double
    a1 = 0.325277,
    a2 = 0.492785,
    a3 = 0.310289,
    a4 = 1.317104,
    a5 = 2.425681;
    
  /* Jenkins et al. SO 180 best-fit
   */
  if(SO180)
    {
      JENKINS_A = 0.301;
      JENKINS_B = 0.64;
      JENKINS_C = 3.82;
    }      
  if(SO324)
    {
      JENKINS_A = 0.316;
      JENKINS_B = 0.67;
      JENKINS_C = 3.82;
    }      
  if(JENKINS) //their .2 fof function
    {
      JENKINS_A = 0.315;
      JENKINS_B = 0.61;
      JENKINS_C = 3.8;
    }      


  /* First normalize the power spectrum
   */
  pnorm=SIGMA_8/sigmac(8.0);
  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm*sigmac(rm);
  logm=log10(mass);
  
  mlo=0.99*mass;
  mhi=1.01*mass;
  rlo=pow(3.0*mlo/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rhi=pow(3.0*mhi/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  slo=pnorm*sigmac(rlo);
  shi=pnorm*sigmac(rhi);
  dsdM=(shi-slo)/(mhi-mlo);

  if(SO324)goto JENKINS_FUNCTION;
  if(SO180)goto JENKINS_FUNCTION;
  if(WARREN)goto WARREN_FUNCTION;
  if(ST)goto ST_FUNCTION;
  if(JENKINS)goto JENKINS_FUNCTION;

  /* Tinker et al. (2008) for SO as a function of DELTA_HALO
   */
  if(DELTA_HALO != prev_delta || prev_cosmo != RESET_COSMOLOGY)
    {
      initialize_mass_function(&a1,&a2,&a3,&a4);
      prev_delta = DELTA_HALO;
      prev_cosmo = RESET_COSMOLOGY;
      //a4*=1.3;
      //a1*=1.15;
      fprintf(stderr,"MF PARAMS for DELTA=%f %f %f %f %f\n",DELTA_HALO,a1,a2,a3,a4);
    }
  
  n = -a1*(pow(sig/a3,-a2)+1)*exp(-a4/sig/sig)*OMEGA_M*RHO_CRIT/mass/sig*dsdM;
  // NB! factor for accounting for satellites
  mass = log10(mass)-11;
  fac = 1.0;
  //fac = 1.1;
  if(mass>0.15)
   fac = pow(mass-0.15,0.5)/(1+exp(pow(mass,1.1))*2)*wpl.satfac + 1;
  fac = 1.0;
  //printf("NB! altering halo mass function by: %f %f %f\n",fac,mass,wpl.satfac);
  n = n*fac;
  return(n);

  /* Jenkins et al. FOF .2 best-fit (unless SO180==1)
   */
 JENKINS_FUNCTION:
  a=-JENKINS_A*OMEGA_M*RHO_CRIT/mass/sig;
  n=a*dsdM*exp(-pow(fabs(JENKINS_B-log(sig)),JENKINS_C));
  return(n);

  /* Warren et al. (calibrated only on concordance cosmology, FOF.2)
   */
 WARREN_FUNCTION:
  n = -0.7234*(pow(sig,-1.625)+0.2538)*exp(-1.198/sig/sig)*OMEGA_M*RHO_CRIT/mass/sig*dsdM;
  return(n);

 ST_FUNCTION:

  /* This is a bunch of Sheth-Tormen stuff.
   */
  nuprime=0.841*DELTA_CRIT/sig;
  nufnu=0.644*(1+1.0/pow(nuprime,0.6))*(sqrt(nuprime*nuprime/2/PI))*exp(-nuprime*nuprime/2);
  //n=RHO_CRIT*OMEGA_M/mass*mass*nufnu*fabs(dsdM);
  n=RHO_CRIT*OMEGA_M/mass*nufnu*fabs(dsdM)/sig;
  return(n);

}


/* It may be a bit costly to run the above function every time you need
 * dn/dM, so here we put the values into an array and then interpolate. 
 *
 * The currentrange of masses calculated is 10^9 to 10^16.7. The tabulation is
 * done in log(M), so the spline interpolation will perform a power-law fit
 * to masses outside this range.
 */
double dndM_interp(double m)
{
  static int flag=0,prev_cosmo=0, n;
  static double *x,*y,*y2;
  int i;
  double dm,max=16.7,min=8,a,m1,m2,dm1;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      n = 200;
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting mass function for %f %f\n",OMEGA_M,SIGMA_8);
      fflush(stdout);

      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dm=(double)(max-min)/n;
      for(i=1;i<=n;++i)
	{
	  x[i]=pow(10.0,min+i*dm);
	  y[i]=log(halo_mass_function(x[i]));
	  //printf("MF%d %e %e\n",RESET_COSMOLOGY,x[i],exp(y[i]));fflush(stdout);
	  if(isnan(y[i])) { n = i-1; break; }
	  if(isinf(y[i])) { n = i-1; break; }
	  x[i]=log(x[i]);
	  continue;
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
      //fprintf(stderr,"MMAX %e\n",exp(x[n]));
    }
  m=log(m);
  if(m>x[n])return 0;
  splint(x,y,y2,n,m,&a);
  return(exp(a));

}


/* Reads mass function in from a file and spline interpolates.
 */
double nbody_massfunction(double m)
{
  static int flag=0,n;
  static double *x,*y,*y2,log10_2,normhi,normlo;
  float x1,x2,x3;
  int i;
  double a,dx;
  char aa[1000];
  FILE *fp;

  if(!flag)
    {
      log10_2=log10(2.0);
      flag++;
      if(!(fp=fopen(Files.MassFuncFile,"r")))
	endrun("ERROR opening MassFuncFile");
      i=0;
      n = filesize(fp);
      fprintf(stderr,"Read %d lines from [%s]\n",n,Files.MassFuncFile);
      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x1,&x2);
	  x[i]=log(x1);
	  y[i]=log(x2);
	  fgets(aa,1000,fp);
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      fclose(fp);
      fprintf(stderr,"Minimum halo mass in N-body dndM= %e\n",exp(x[1]));

      normhi = exp(y[n])/halo_mass_function(exp(x[n]));
      normlo = exp(y[1])/halo_mass_function(exp(x[1]));
    }

  m=log(m);
  //  if(m>x[n])return(0);
  /*if(m<x[1])return(0);*/

  splint(x,y,y2,n,m,&a);
  return(exp(a));
}      

void initialize_mass_function(double *a1, double *a2, double *a3, double *a4)
{
  int n = 9, i;
  double *x, *y, *z, at, ztemp;

  x = dvector(1,n);
  y = dvector(1,n);
  z = dvector(1,n);

  // initialize the overdensities
  for(i=1;i<=9;i+=2)
    x[i] = log(200*pow(2.0,(i-1.0)/2.0));
  for(i=2;i<=9;i+=2)
    x[i] = log(300*pow(2.0,(i-2.0)/2.0));

  //first parameter
  y[1] = 1.858659e-01 ;
  y[2] = 1.995973e-01 ;
  y[3] = 2.115659e-01 ;
  y[4] = 2.184113e-01 ;
  y[5] = 2.480968e-01 ;
  y[6] = 2.546053e-01 ;
  y[7] = 2.600000e-01 ;
  y[8] = 2.600000e-01 ;
  y[9] = 2.600000e-01 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a1);
  if(DELTA_HALO>=1600) *a1 = 0.26;

  //second parameter
  y[1] = 1.466904e+00 ;
  y[2] = 1.521782e+00 ;
  y[3] = 1.559186e+00 ;
  y[4] = 1.614585e+00 ;
  y[5] = 1.869936e+00 ;
  y[6] = 2.128056e+00 ;
  y[7] = 2.301275e+00 ;
  y[8] = 2.529241e+00 ;
  y[9] = 2.661983e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a2);

  //third parameter
  y[1] = 2.571104e+00 ;
  y[2] = 2.254217e+00 ;
  y[3] = 2.048674e+00 ;
  y[4] = 1.869559e+00 ;
  y[5] = 1.588649e+00 ;
  y[6] = 1.507134e+00 ;
  y[7] = 1.464374e+00 ;
  y[8] = 1.436827e+00 ;
  y[9] = 1.405210e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a3);


  //fourth parameter
  y[1] = 1.193958e+00;
  y[2] = 1.270316e+00;
  y[3] = 1.335191e+00;
  y[4] = 1.446266e+00;
  y[5] = 1.581345e+00;
  y[6] = 1.795050e+00;
  y[7] = 1.965613e+00;
  y[8] = 2.237466e+00;
  y[9] = 2.439729e+00;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a4);

  // now adjust for redshift
  if(!(REDSHIFT>0))return;

  ztemp = REDSHIFT;
  if(REDSHIFT>3) ztemp = 3.0;
  *a1 *= pow(1+ztemp,-0.14);
  *a2 *= pow(1+ztemp,-0.14);
  at = -pow(0.75/log10(DELTA_HALO/75),1.2);
  at = pow(10.0,at);
  *a3 *= pow(1+ztemp,-at);

}

/* This is using the normalized mass functions that
 * integrate to the mean matter density of the universe
 */
void initialize_mass_function2(double *a1, double *a2, double *a3, double *a4, double *a5)
{
  int n = 9, i;
  double *x, *y, *z;

  x = dvector(1,n);
  y = dvector(1,n);
  z = dvector(1,n);

  // initialize the overdensities
  for(i=1;i<=9;i+=2)
    x[i] = log(200*pow(2.0,i-1.0));
  for(i=2;i<=9;i+=2)
    x[i] = log(300*pow(2.0,i-1.0));

  //first parameter
  y[1] = 5.143873e-01 ;
  y[2] = 4.782773e-01 ;
  y[3] = 4.773988e-01 ;
  y[4] = 4.548802e-01 ;
  y[5] = 4.431224e-01 ;
  y[6] = 3.982487e-01 ;
  y[7] = 3.976351e-01 ;
  y[8] = 3.616926e-01 ;
  y[9] = 3.282215e-01 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a1);

  //second parameter
  y[1] = 1.973173e+00 ;
  y[2] = 2.056993e+00 ;
  y[3] = 2.297465e+00 ;
  y[4] = 2.563651e+00 ;
  y[5] = 2.834270e+00 ;
  y[6] = 2.916583e+00 ;
  y[7] = 3.294007e+00 ;
  y[8] = 3.369249e+00 ;
  y[9] = 3.300715e+00 ;
  
  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a2);

  //third parameter
  y[1] = 9.951248e-01 ;
  y[2] = 9.895878e-01 ;
  y[3] = 9.343687e-01 ;
  y[4] = 9.303064e-01 ;
  y[5] = 9.588204e-01 ;
  y[6] = 1.044357e+00 ;
  y[7] = 1.065549e+00 ;
  y[8] = 1.119347e+00 ;
  y[9] = 1.164469e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a3);

  //fourth parameter
  y[1] = 1.227880e+00 ;
  y[2] = 1.310414e+00 ;
  y[3] = 1.403302e+00 ;
  y[4] = 1.553185e+00 ;
  y[5] = 1.701724e+00 ;
  y[6] = 1.906987e+00 ;
  y[7] = 2.137733e+00 ;
  y[8] = 2.394142e+00 ;
  y[9] = 2.572204e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a4);


  //normalization parameter
  y[1] = 4.816036e-01 ;
  y[2] = 4.660390e-01 ;
  y[3] = 4.935099e-01 ;
  y[4] = 4.937660e-01 ;
  y[5] = 4.960265e-01 ;
  y[6] = 4.496969e-01 ;
  y[7] = 4.663272e-01 ;
  y[8] = 4.288530e-01 ;
  y[9] = 3.879745e-01 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a5);

  //set up mass funciton parameters as globals
  /*
  MFPARAM1 = *a1;
  MFPARAM2 = *a2;
  MFPARAM3 = *a3;
  MFPARAM4 = *a4;
  MFPARAM5 = *a5; 
  */
}
