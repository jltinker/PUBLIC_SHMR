#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

double rs_g7, rs2_g7, rp_g7, n_g7, ncen_g7, nsat_g7, cvir_g7, rvir_g7, n_g7,
  cvir_g7, rvir_g7, sigv_g7, sigvc_g7, moment_ss_g7;

double func1_zspace(double m);
double func2_zspace(double m);


double one_halo(double rs, double rp)
{
  int i;
  double s1,rlim,zlim,t0,t1,r,mlo,mhi,dm,ma,mb;

  

  if(!HOD.pdfs)return(0);

  if(OUTPUT>1) {
    printf("xi1h_zspace> rs= %.2f rp= %.2f\n",rs,rp);
    fflush(stdout);
  }

  rs_g7=rs;
  rs2_g7 = rs*rs;
  rp_g7=rp;
  r=sqrt(rs*rs+rp*rp);

  /* The smallest halo allowed in the integral is where r_vir(M) = rs/2
   */
  mlo = 4./3.*DELTA_HALO*PI*rs*rs*rs*OMEGA_M*RHO_CRIT*0.125;
  if(mlo<HOD.M_low)mlo=HOD.M_low;
  if(mlo<HOD.M_cut && (HOD.pdfs==2 || HOD.pdfs==5))mlo = HOD.M_cut;

  mlo=log(mlo);
  mhi=log(HOD.M_max);

  /* If the velocity bias is zero, then the one-halo term
   * is the same in real and redshift space.
   */
  if(VBIAS==0) { return(one_halo_real_space(r)); muh(999); exit(0); }

  /* Check to make sure we're not totally out of the range of the one-halo term.
   */
  rlim=2*pow(3*HOD.M_max/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),THIRD);
  if(rlim<rs)
    return(0);

  t0 = second();
  s1=qromo(func1_zspace,mlo,mhi,midpnt)/(2*PI*GALAXY_DENSITY*GALAXY_DENSITY);
  // The below 2 lines are for fiber collisions
  //if(rs<0.1)s1*=0.15;
  //if(rs<0.2 && rs>=0.1)s1*=0.75;
  t1 = second();
  if(OUTPUT>1) {
    printf("xi1h_zspace> TIME %.2f %e\n",timediff(t0,t1),s1);
    printf("xi1h_zspace> rs= %.2f rp= %.2f xi= %e\n",rs,rp,s1);
    fflush(stdout);
  }
  return(s1);

}

/* This is the integral over dM which calls the integral over dz
 */
double func1_zspace(double m)
{
  static int prev_cosmo=0;
  double s1,s2,zlim,rlim,logm,vbias,z,t0,t1;
  static double fac=-1;
  int i;

  logm=m;
  m=exp(m);
  rlim=pow(3*m/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),THIRD);
  rvir_g7=2*rlim;
  rlim*=2;
  if(rs_g7>rlim)return(0);
  zlim=sqrt(rlim*rlim-rs_g7*rs_g7);

  /* Put a few things in globals
   */
  cvir_g7 = halo_concentration(m)*CVIR_FAC;
  n_g7 = dndM_interp(m);
  nsat_g7 = N_sat(m);
  ncen_g7 = N_cen(m);
  moment_ss_g7=moment_ss(m);

  vbias=VBIAS;

  /* TEMP-> Testing a vbias at some mass threshold (or with some mass dependence)
   */
  /*
  vbias = (VBIAS_SLOPE*(logm-34.5)+VBIAS);
  if(vbias<0.1)vbias=0.1;

  vbias=1;
  if(VBIAS_MASS_THRESHOLD>0)
    if(m>VBIAS_MASS_THRESHOLD)
      vbias=VBIAS;
  */
  //  vbias = 0.4/3*logm/LOGE_10 - 1;
  //if(vbias<0.6)vbias=0.6;

  if(fac<0 || RESET_COSMOLOGY!=prev_cosmo)
    {
      fac=sqrt(4.499E-48/2.0)*
	pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3.,1.0/6.0)*3.09E19;
      fac *= sqrt(1+REDSHIFT);
    }

  sigv_g7=vbias*fac*pow(m,0.33333333333);

  if(VBIAS_C)
    sigvc_g7=(fac*pow(m,0.33333333333))*sqrt(VBIAS_C*VBIAS_C + vbias*vbias);
  else
    sigvc_g7=sigv_g7;

  // printf("SIGV %e %e %e\n",m,sigv_g7,sigvc_g7);

  /*
  t0 = second();
  for(i=1;i<=100;++i)
    {
      z = (2*zlim)/100*(i-0.5) - zlim;
      s2 += func2_zspace(z)*(zlim/50.)*100;
    }
  t1 = second();
  printf("BOO %.2f %e %e\n",timediff(t0,t1),m,s2*m);
  return(s2*m);
  */
  
  s2=qtrap(func2_zspace,-zlim,zlim,1.0E-4)*100*HUBBLEZ; //-6 solved for Fisher
  /* s2=qromo(func2_zspace,-zlim,zlim,midpnt)*100;  */
  prev_cosmo=RESET_COSMOLOGY;
  return(s2*m);
}

double func2_zspace(double z)
{
  double N,n,rvir,f_ss,f_cs=0,cvir,x,sigv,sigvc,vz,r;

  r=sqrt(rs2_g7+z*z);
  rvir=rvir_g7;
  cvir=cvir_g7;
  sigv=sigv_g7;
  sigvc=sigvc_g7;
  n=n_g7;

  vz=100*(rp_g7-z)*HUBBLEZ;

  f_ss=dFdx_ss_interp(r/rvir,cvir)*moment_ss_g7*
    0.5/RT2PI*exp(-vz*vz/(4.0*sigv*sigv))/(sigv*ROOT2);
  if(rs_g7<=rvir/2)
    f_cs=dFdx_cs(r/rvir,cvir)*nsat_g7*ncen_g7/
      RT2PI*exp(-vz*vz/(2.0*sigvc*sigvc))/sigvc;

  x=n*(f_ss+f_cs)/(rvir*r*r);
  return(x);

}

