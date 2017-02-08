#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"

/* Local variables for the coordinate in redshift space.
 */
double rs_g8,rp_g8,
  iorder_k1, r_k1;

/* Internal functions.
 */
double func_kaiser(double z);
double xi_bar(double r);
double f_xibar(double r);
double xi_prime(double r, double mu);
double xi_t(double r);
double Lpoly(double x, int i);
double f_xi2bar(double r);
double xi_harmonic(double r, int i);
double scale_dependent_dispersion(double r);
double scale_dependent_skewness(double r);
double func_xip(double k);
double xi_p(double r, int iorder);
double func1k(double k);
double xi_harmonic_pk(double r, int i);

void initial_dispersion_model_values(double *a, double **pp, double *yy);
void dispersion_model_input(void);
double chi2_dispersion_model(double *a);

/* Internal variables.
 */
struct DISPERION_MODEL_STRUCT {
  int n;
  double **rs,**rp;
  double **x,**e;
} XX;

double  *rr_g8, *aa_g8;
int nparams_g8, niter_g8=0;

void fit_dispersion_model()
{
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m;
  FILE *fp;
  char aa[1000];

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF xi(SIGMA,PI) WITH DISPERSION MODEL.\n");
  fprintf(stderr,    "--------------------------------------------------------\n\n");

  if(POWELL)
    FTOL=1.0E-4;
  else
    FTOL=1.0E-4;

  nparams_g8 = n = 8;

  dispersion_model_input();

  wp.ncf=n;
  a=dvector(1,n);
  rr_g8 = dvector(1,n);
  aa_g8 = dvector(1,n);

  for(i=1;i<=n;++i)
    rr_g8[i] = (log(10)-log(0.1))/(n-1.)*(i-1) + log(0.1);

  BETA = 0.5;

  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  initial_dispersion_model_values(a,pp,yy);

  if(POWELL) 
    {
      if(OUTPUT)printf("dispersion_model> starting powell.\n");
      powell(a,pp,n,FTOL,&niter,&chi2min,chi2_dispersion_model);
      chi2min = chi2_dispersion_model(a);
    }
  else
    {
      if(OUTPUT)printf("dispersion_model> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_dispersion_model,&niter);
      for(i=1;i<=n;++i)a[i]=pp[1][i];
      chi2min = chi2_dispersion_model(a);
    }	
}


void initial_dispersion_model_values(double *a, double **pp, double *yy)
{
  static int flag=0;
  int i,j;
  double d[100];


  if(!flag) {
    for(i=1;i<=nparams_g8;++i)
      a[i] = 500;

    printf("INITIAL VALUES: ");
    for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
    printf("\n");
  }
  flag++;

  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=wp.ncf;++i)
    d[i]=a[i]*0.1/flag;


  if(POWELL)
    {
      for(i=1;i<=wp.ncf;++i)
	{
	  for(j=1;j<=wp.ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(j=1;j<=wp.ncf;++j)
	pp[1][j]=a[j];
      yy[1]=chi2_dispersion_model(a);
    
      for(i=1;i<=wp.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	    yy[i+1]=chi2_dispersion_model(a);	  
	  for(j=1;j<=wp.ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[wp.ncf]-=d[wp.ncf];
    }
}

void dispersion_model_input()
{
  FILE *fp;
  int i,j;

  fp = openfile("xi2d.data");
  XX.n = (int)sqrt(filesize(fp)*0.99)+1;
  
  XX.rs = dmatrix(1,XX.n,1,XX.n);
  XX.rp = dmatrix(1,XX.n,1,XX.n);
  XX.x  = dmatrix(1,XX.n,1,XX.n);
  XX.e  = dmatrix(1,XX.n,1,XX.n);

  for(i=1;i<=XX.n;++i)
    for(j=1;j<=XX.n;++j)
      fscanf(fp,"%lf %lf %lf %lf",&XX.rs[i][j],&XX.rp[i][j],&XX.x[i][j],&XX.e[i][j]);

  fclose(fp);

}

double chi2_dispersion_model(double *a)
{
  int i,j;
  double xx,r,chi2=0;

  for(i=1;i<=nparams_g8;++i)
    {
      aa_g8[i] = a[i];
      if(a[i]<=0)return(1.0E7);
    }

  for(i=2;i<=XX.n;++i)
    for(j=2;j<=XX.n;++j)
      {
	xx = kaiser_distortion(XX.rs[i][j],XX.rp[i][j]);
	/* printf("%f %f %f %f %f\n",XX.rs[i][j],XX.rp[i][j],XX.x[i][j],xx,xx); */
	chi2 += (xx-XX.x[i][j])*(xx-XX.x[i][j])/(XX.e[i][j]*XX.e[i][j]);
      }
  niter_g8++;
  printf("ITER %d %e ",niter_g8,chi2);
  for(i=1;i<=nparams_g8;++i)
    printf("%e ",a[i]);
  printf("\n");
  fflush(stdout);

  if(niter_g8%100==0)
    {
      for(i=1;i<=100;++i)
	{
	  r = exp((log(40)-log(0.1))/(100-1.)*(i-1) + log(0.1));
	  xx = scale_dependent_dispersion(r);
	  printf("PVD %d %f %f\n",niter_g8,r,xx);
	}
      fflush(stdout);
    }

  return(chi2);
}



/*************************************************************
 * This is the function called in place of the full 2HALO term
 *************************************************************/

double kaiser_distortion(double rs, double rp)
{
  double s1,rlim=100,rr;
  static int prev=-1;

  rs_g8=rs;
  rp_g8=rp;
  rr=sqrt(rs*rs+rp*rp);
  rlim = 20*SIGV/(100.0*HUBBLEZ); // go out 4 sigma around location
  s1=qtrap(func_kaiser,-rlim,rlim,1.0E-5); 
  return(s1-1);
}

double linear_kaiser_distortion(double r, double z)
{
  double mu, x;

  mu= sqrt(1 - z/r*z/r);
  mu = z/r;
  x = xi_prime(r,mu);
  if(x<-1)return -1;
  return x;
}

double func_kaiser(double z)
{
  double r,mu,v,x,sigv;

  // new version
  r=sqrt((z+rp_g8)*(z+rp_g8)+rs_g8*rs_g8);
  mu=(z+rp_g8)/r;
  v=z*100.0*HUBBLEZ;

  sigv = SIGV*pow(OMEGA_M/0.3,GAMMA)*(SIGMA_8/0.76);
  // NB! tryingt to shift the dispersion by scale
  //sigv = scale_dependent_dispersion(r)*pow(OMEGA_M/0.3,GAMMA)*(SIGMA_8/0.76);
  x=(1+xi_prime(r,mu))*exp(-ROOT2*fabs(v)/sigv)/ROOT2/sigv*100.0*HUBBLEZ;
  return x;


  // original version
  r=sqrt(z*z+rs_g8*rs_g8);
  mu=z/r;
  v=(rp_g8-z)*100.0*HUBBLEZ;

  sigv = SIGV;
  x=(1+xi_prime(r,mu))*exp(-ROOT2*fabs(v)/sigv)/ROOT2/sigv*100.0*HUBBLEZ;
  return x;




  sigv = scale_dependent_dispersion(r);
  /*
  if(v*z<0)
    sigv*=scale_dependent_skewness(r);
  */
  /*
  x=(1+one_halo_real_space(r)+two_halo_real_space(r))*
    exp(-ROOT2*fabs(v)/sigv)/ROOT2/sigv*100.0;
  */
  x=(1+xi_prime(r,mu)+one_halo_real_space(r))*exp(-ROOT2*fabs(v)/sigv)/ROOT2/sigv*100.0*HUBBLEZ;
  return(x);
}
double scale_dependent_skewness(double r)
{
  return 1.5;
}

double scale_dependent_dispersion(double r)
{
  static int niter=-1,flag=0;
  static double *z,*x,*y;
  double a;
  int i;

  return SIGV*(-0.5*log(r/10+1)+(1+0.5*log(2)));

  if(niter!=niter_g8)
    {
      niter = niter_g8;
      if(!flag)
	{
	  flag++;
	  x=dvector(1,nparams_g8+1);
	  y=dvector(1,nparams_g8+1);
	  z=dvector(1,nparams_g8+1);
	  y[1] = 0;
	  x[1] = log(0.01);
	  for(i=2;i<=nparams_g8+1;++i)
	    x[i] = rr_g8[i-1];
	}
      for(i=2;i<=nparams_g8+1;++i)
	y[i] = aa_g8[i-1];
      spline(x,y,nparams_g8+1,1.0E+30,1.0E+30,z);
    }
  r=log(r);
  if(r>rr_g8[nparams_g8])return(aa_g8[nparams_g8]);
  splint(x,y,z,nparams_g8+1,r,&a);
  return(a);
}

/* This function tabulates the spherically averaged correlation function.
 * Equation (A7) in Madgwick et al.
 */
double xi_bar(double r)
{
  static int flag=0,prev=-1;
  static double *x,*y,*y2;
  int n=100,i;
  static double rmin,rmax,lgr,a;

  if(!flag || prev!=RESET_KAISER)
    {
      if(OUTPUT>1)
	printf("Tabulating xi_bar...\n");
      prev=RESET_KAISER;
      if(!flag)
	{
	  x=malloc(n*sizeof(double));
	  x--;
	  y=malloc(n*sizeof(double));
	  y--;
	  y2=malloc(n*sizeof(double));
	  y2--;
	}
      flag=1;
      rmin=-1.9;
      if(R_MIN_2HALO>0)
	rmin=log10(R_MIN_2HALO);
      rmax=log10(R_MAX_2HALO);
      for(i=1;i<=n;++i)
	{
	  lgr=(rmax-rmin)*(i-1.0)/(n-1.0)+rmin;
	  x[i]=pow(10.0,lgr);
	  y[i]=qtrap(f_xibar,0.0,x[i],1.0E-4)*3/(x[i]*x[i]*x[i]);
	  //printf("XIBAR %f %e %f\n",x[i],y[i]);
	}
      spline(x,y,n,1.0E30,1.0E30,y2);
    }
  
  if(r<x[1])return(-1);
  splint(x,y,y2,n,r,&a);
  //if(a<xi_t(r))return(xi_t(r));
  return(a);
}

double f_xibar(double r)
{
  double x;
  x = xi_t(r);
  //x=two_halo_real_space(r);
  /* x=one_halo_real_space(r)+two_halo_real_space(r); */
  return(x*r*r);
}

/* This function tabulates the xi(r)*r^4 dr function.
 * Equation (A8)
 */
double xi_2bar(double r)
{
  static int flag=0,prev=-1;
  static double *x,*y,*y2,rmin;
  int n=100,i;
  double rmax,lgr,a;

  if(!flag || prev!=RESET_KAISER)
    {
      if(OUTPUT>1)
	printf("Tabulating xi_2bar...\n");
      prev=RESET_KAISER;
      if(!flag)
	{
	  x=malloc(n*sizeof(double));
	  x--;
	  y=malloc(n*sizeof(double));
	  y--;
	  y2=malloc(n*sizeof(double));
	  y2--;
	}
      flag=1;
      rmin=-1.9;
      if(R_MIN_2HALO>0)
	rmin = log10(R_MIN_2HALO);
      rmax=log10(R_MAX_2HALO);
      for(i=1;i<=n;++i)
	{
	  lgr=(rmax-rmin)*(i-1.0)/(n-1.0)+rmin;
	  x[i]=pow(10.0,lgr);
	  y[i]=qtrap(f_xi2bar,0.0,x[i],1.0E-4)*5/(x[i]*x[i]*x[i]*x[i]*x[i]);
	  //printf("XI2BAR %f %e\n",x[i],y[i]);
	}
      spline(x,y,n,1.0E30,1.0E30,y2);
    }

  if(r<x[1])return(-1);
  splint(x,y,y2,n,r,&a);
  //  if(a<xi_t(r))return(xi_t(r));
  return(a);
}

double f_xi2bar(double r)
{
  double x;
  x = xi_t(r);
  //x=two_halo_real_space(r);
  /* x=one_halo_real_space(r)+two_halo_real_space(r); */
  return(x*r*r*r*r);
}

/* These are Legendre polynomials
 */
double Lpoly(double x, int i)
{
  switch(i) {
  case 0: return(1.0);
  case 2: return(0.5*(3*x*x-1));
  case 4: return((35*x*x*x*x-30*x*x+3)/8.0);
    /*case 4: return((384*x*x*x*x+1152*x*x*(x*x-1)+144*(x*x-1)*(x*x-1))/(16.*24.));*/
  default: return(0.0);
  }
}

double xi_harmonic(double r, int i)
{
  switch(i){
  case 0: return((1+2.*BETA/3.+BETA*BETA/5.)*xi_t(r));
  case 2: return((4.*BETA/3.+4.*BETA*BETA/7.)*(xi_t(r)-xi_bar(r)));
  case 4: return(8.*BETA*BETA/35.*(xi_t(r)+2.5*xi_bar(r)-3.5*xi_2bar(r)));
  default: return(0);
  }
}

double xi_prime(double r, double mu)
{
  double x0,x2,x4;
  
  x0=xi_harmonic(r,0)*Lpoly(mu,0);
  x2=xi_harmonic(r,2)*Lpoly(mu,2);
  x4=xi_harmonic(r,4)*Lpoly(mu,4);
  /*
  x0=xi_harmonic_pk(r,0)*Lpoly(mu,0);
  x2=xi_harmonic_pk(r,2)*Lpoly(mu,2);
  x4=xi_harmonic_pk(r,4)*Lpoly(mu,4);
  */
  return(x0+x2+x4);
}

double xi_t(double r)
{
  return(two_halo_real_space(r));
  return(one_halo_real_space(r)+two_halo_real_space(r));
}


