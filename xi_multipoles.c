
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* This is a simple integrator to calculate the correlation function multipoles.
 * The xi(sigma,pi) will be calculated on a polar grid in log r, phi, and then
 * simple midpoint method to integrate over angle at each r.
 * There will be 16 angular bins and maybe 20 radial bins.
 *
 * I suppose I could interpolate over angle from the 16 calculations, but this is the
 * most direct way of comparing to the n-body simulations, which are calculated exactly
 * this way.
 */

int nrad_g6,nphi_g6,NR,nrad_g6_temp=50;
double *phi_g6,*rr_g6,*xi_mono,**xi_g6,*xi_mono_temp,*rr_g6_temp;
double spherically_averaged_xi(double s);
double spherically_averaged_xi_temp(double s);
double monopole(double mu);
double quadrupole(double mu);

void xi_multipoles(double rr, double *xx1, double *xx2)
{
  static int flag=1, niter=1;
  static double *xi_quad, *xi_quad_zz, *xi_mono_zz;

  int i,j,k,n,istep=1,istart=1,jstep=1,jstart=1,sdss_ilimit=3;
  double rlo,rhi,rs,rp,delta_phi,phi1,t1h=0,t2h=0,xi,
    x1,x2,x3,x1h,x2h,a,b,t0,t1,*send,*recv,dlogr,xi_bar,drbin,
    philo,phihi,rp_min=0;
  FILE *fp;
  char fname[100],aa[100];
  static int flag1 = 1;

#ifdef PARALLEL
  istep = NTask;
  istart = ThisTask + 1;
  jstart = 1;
  jstep = 1;
#endif

  muh(RESET_ZSPACE);
  fmuh(SIGV);
  muh(KAISER);
  if(!RESET_ZSPACE) goto INTERP_ONLY;


  // linear spacing
  nrad_g6=40;
  rlo = 0;
  rhi = 40.0;
  dlogr = 1; 

  // log spacing
  nrad_g6=25;
  rlo = 0.1;
  rhi = 40.0;
  dlogr = log(rhi/rlo)/(nrad_g6-1);



  nphi_g6 = 18;
  delta_phi = PI/2./nphi_g6;

  if(flag) {
    flag=0;
    rr_g6=dvector(1,nrad_g6);
    phi_g6=dvector(1,nphi_g6);
    xi_mono=dvector(1,nrad_g6);
    xi_quad=dvector(1,nrad_g6);
    xi_g6=dmatrix(1,nrad_g6,1,nphi_g6);
  }

  
  recv=dvector(1,nrad_g6);
  for(i=1;i<=nrad_g6;++i)
    xi_mono[i]=xi_quad[i]=0;

  x1=one_halo_real_space(1)+two_halo_real_space(1);
  //muh(999);
  x1 = two_halo(50,50);
  //exit(0);

  t0=second();

  //logbins
  rlo = 0.1;
  dlogr = log(40/0.1)/(nrad_g6-1);

  /*
  for(i=-10;i<=10;++i)
    {
      x1 = pow(10.0,i/10.0);
      printf("%e %e %e\n",x1,one_halo(0.1,x1),one_halo_real_space(sqrt(x1*x1+0.1*0.1)));
    }
  exit(0);
  */	

  //sprintf(fname,"multipoles.%d",niter++);
  sprintf(fname,"%s.multipoles",Task.root_filename);
  fp = fopen(fname,"w");

  /* Calculate the redshift space correlation function on the polar grid.
   */
  for(i=istart;i<=nrad_g6;i+=istep)
    {

      NR=i;
      xi_mono[i]=0;
      xi_quad[i]=0;

      rr_g6[i] = exp((i-1)*dlogr)*rlo;
      //rr_g6[i] = (i-0.5)*dlogr + rlo; // linear bins

      for(j=nphi_g6;j>=1;--j)
	{
	  phi_g6[j]=(j-0.5)/(nphi_g6)*PI/2.0;
	  rs=rr_g6[i]*cos(phi_g6[j]);
	  rp=rr_g6[i]*sin(phi_g6[j]);

	  t0=second();
	  x1h=one_halo(rs,rp);
	  t1=second();
	  t1h+=timediff(t0,t1);
	  
	  t0=second();
	  x2h=two_halo(rs,rp);
	  t1=second();
	  t2h+=timediff(t0,t1);
	  
	  xi=x1h+x2h;
	  xi_g6[i][j] = xi;
	  
	  t1=second();
	  
	  /* Calculate the monopole and quadrupole for this bin.
	   */
	  phi1=PI/2 - phi_g6[j];

	  xi_mono[i]+=xi*sin(phi1)*delta_phi;
	  xi_quad[i]+=5*(3*cos(phi1)*cos(phi1)-1)/2.0*sin(phi1)*delta_phi*xi; 
	  // printf("%f %e %f %f %e %e\n",phi_g6[j],xi,rs,rp,x1h,sin(phi1));
	}
      //exit(0);

      /* Use spline interpolation to integrate the multipoles
       * unless in the small-r regime.
       */
      for(j=1;j<=nphi_g6;++j)
	phi_g6[j]=cos(PI/2 - phi_g6[j]);

      for(xi_bar=0,j=1;j<=i;++j)
	xi_bar+=xi_mono[j]*rr_g6[j]*rr_g6[j]*rr_g6[j]*dlogr;
      xi_bar*=(3/pow(rr_g6[i],3.0));

      if(OUTPUT>=0)
	{
	  fprintf(fp,"MULTIPOLES%d %f %e %e %e %.2f %.2f %.3f %.3f %.3f\n",
		  ThisTask,rr_g6[i],xi_mono[i],xi_quad[i],-xi_quad[i]/(xi_bar - xi_mono[i]),t1h,t2h,
		  xi_mono[i]/(one_halo_real_space(rr_g6[i]) + two_halo_real_space(rr_g6[i])),1+2./3.*BETA+0.2*BETA*BETA,
		  (4./3.*BETA+4./7*BETA*BETA)/(1+2./3.*BETA+0.2*BETA*BETA));
	  fflush(fp);
	  fprintf(stdout,"MULTIPOLES%d %f %e %e %e %.2f %.2f %.3f %.3f %.3f\n",
		  ThisTask,rr_g6[i],xi_mono[i],xi_quad[i],-xi_quad[i]/(xi_bar - xi_mono[i]),t1h,t2h,
		  xi_mono[i]/(one_halo_real_space(rr_g6[i]) + two_halo_real_space(rr_g6[i])),1+2./3.*BETA+0.2*BETA*BETA,
		  (4./3.*BETA+4./7*BETA*BETA)/(1+2./3.*BETA+0.2*BETA*BETA));
	  fflush(stdout);
	}
      // replace with Q(r)
      xi_quad[i] = xi_quad[i]/(xi_mono[i]-xi_bar);

    }
  fclose(fp);
  /*
  // get the quadrupole moment
  for(i=1;i<=nrad_g6;++i)
    {
      xi_bar = spherically_averaged_xi(rr_g6[i])*(3/pow(rr_g6[i],3.0));
      xi_quad[i] = -xi_quad[i]/(xi_bar - xi_mono[i]);
      printf("%e %e %e %e\n",rr_g6[i],xi_mono[i],xi_bar,xi_quad[i]);
    }
  */

  // set up interpolation arrays
  xi_mono_zz = dvector(1,nrad_g6);
  xi_quad_zz = dvector(1,nrad_g6);
  spline(rr_g6, xi_mono, nrad_g6, 1.0E+30, 1.0E+30,xi_mono_zz);
  spline(rr_g6, xi_quad, nrad_g6, 1.0E+30, 1.0E+30,xi_quad_zz);
  
 INTERP_ONLY:
  if(RESET_ZSPACE)
    RESET_ZSPACE = 0 ;

  splint(rr_g6, xi_mono, xi_mono_zz, nrad_g6, rr, xx1);
  splint(rr_g6, xi_quad, xi_quad_zz, nrad_g6, rr, xx2);

}

/* Use the tabluated values of xi(s,p) and calculate the monpole.
 */
double monopole(double mu)
{
  static double *y2;
  static int flag=0, prev_nr=0;
  double a;

  if(prev_nr!=NR)
    {
      y2=dvector(1,nphi_g6);
      spline(phi_g6,xi_g6[NR],nphi_g6,1.0E+30,1.0E+30,y2);
    }
  prev_nr=NR;
  splint(phi_g6,xi_g6[NR],y2,nphi_g6,mu,&a);
  return(a);
}

/* Use the tabluated values of xi(s,p) and calculate the quadrupole.
 */
 
double quadrupole(double mu)
{
  static double *y2;
  static int flag=0, prev_nr=0;
  double a;
  if(prev_nr!=NR)
    {
      y2=dvector(1,nphi_g6);
      spline(phi_g6,xi_g6[NR],nphi_g6,1.0E+30,1.0E+30,y2);
    }
  prev_nr=NR;
  splint(phi_g6,xi_g6[NR],y2,nphi_g6,mu,&a);
  return(2.5*(3*mu*mu-1)*a);
}


/* For the quadrupole to monopole ratio, calculate the spherically averaged
 * monopole as in Hawkins etal. (2003)
 * For Kaiser model, FLAG=0. For HOD model, FLAG=1
 * LINEAR SPACING IN R!!
 */
double spherically_averaged_xi(double s)
{
  static double *y2,*k2;
  static int flag=0,prev=0;
  double a;

  if(!flag || prev!=RESET_COSMOLOGY)
    {
      if(!flag) {
	y2=dvector(1,nrad_g6);
	k2=dvector(1,nrad_g6); 
      }
      flag++;
      spline(rr_g6,xi_mono,nrad_g6,1.0E+30,1.0E+30,y2);
      prev=RESET_COSMOLOGY;
    }
  splint(rr_g6,xi_mono,y2,nrad_g6,s,&a);    
  return(s*s*a);
}

double spherically_averaged_xi_temp(double s)
{
  static double *y2,*k2;
  static int flag=0,prev=0;
  double a;

  if(!flag || prev!=RESET_COSMOLOGY)
    {
      if(!flag) {
	y2=dvector(1,nrad_g6_temp);
	k2=dvector(1,nrad_g6_temp); 
      }
      flag++;
      spline(rr_g6_temp,xi_mono_temp,nrad_g6_temp,1.0E+30,1.0E+30,y2);
      prev=RESET_COSMOLOGY;
    }
  splint(rr_g6_temp,xi_mono_temp,y2,nrad_g6_temp,s,&a);    
  return(s*s*a);
}









