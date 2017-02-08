#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

double rs_g5, rs2_g5, rp_g5;

double streaming_model_integrand(double z);
int xi2h_zero_check(double z1, double z2);


/* The redshift-space correlation function evalated at r_sig, r_pi:
 *
 * This calculation is done through direct implementation of the
 * streaming model on the two-halo term:
 *
 *    xi_2h(rs,rp) = \int [1+xi_2h(r)] P(vz) dz
 *
 * where: r^2 = rs^2 + z^2
 *        vz = 100*(rp - z)
 *        P(vz) = the N_gal(M) weighted halo pairwise velocity distribution,
 *              tabulated elsewhere.
 */
double two_halo(double rs, double rp)
{
  double r,zmax,zmin,s1,z2h_min,z2h_max;
  int i;

  if(KAISER)
    return kaiser_distortion(rs,rp);

  // NB!! RETURNING linear+exp model 
  //return linear_hod_two_halo_zspace(rs,rp);

  rs_g5=rs;
  rs2_g5=rs*rs;
  rp_g5=rp;

  if(XI_MAX_RADIUS<=0)two_halo_real_space(1.0);

  // Just trying this out for speed purposes
  //  return(linear_kaiser_distortion(sqrt(rs*rs + rp*rp),rp));


  if(OUTPUT>1) {
    printf("xi2h_zspace> rs= %.2f rp= %.2f\n",rs,rp);
    fflush(stdout);
  }

  /* The maximum velocity in the table is ~4000 km/s, so limit
   * the integration to that range in z.
   */
  zmin = rp - mabs(MAXVEL)*0.01;
  zmax = rp + mabs(MAXVEL)*0.01;


  /* Then check to make sure that we're not past the limit
   * of the read-space calculation.
   */
  z2h_max = sqrt(XI_MAX_RADIUS*XI_MAX_RADIUS - rs*rs)*0.99;
  if(zmax>z2h_max)zmax = z2h_max;
  if(zmin<-z2h_max)zmin = -z2h_max;

  /*
  for(i=-200;i<=200;++i)
    {
      rs =  4.290023e-01;
      rp = 0.1652;
      rs_g5=rs;
      rs2_g5=rs*rs;
      rp_g5=rp;
      
      printf("BLAH %f %e\n",i/10.0,streaming_model_integrand(i/10.0));
    }
  exit(0);
  */

  /* If outside 2halo exclusion region, do contiguous integration.
   */
  zmin = -60;
  zmax = 60;
  if(rs>R_MIN_2HALO)
    {
      //muh(100);
      s1 = qromo(streaming_model_integrand,zmin,zmax,midpnt);
      return(s1*100*HUBBLEZ-1);
    }

  /* If inside exclusion region, break up integral in 2 parts.
   */
  z2h_min = sqrt(R_MIN_2HALO*R_MIN_2HALO - rs*rs);
  if(zmin<z2h_min)
    s1 = qromo(streaming_model_integrand,z2h_min,zmax,midpnt);
  else
    s1 = qromo(streaming_model_integrand,zmin,zmax,midpnt);
  
  if(zmin<-z2h_min)
    if(xi2h_zero_check(zmin,-z2h_min))
      s1+= qromo(streaming_model_integrand,zmin,-z2h_min,midpnt);
  
  return(s1*100*HUBBLEZ-1);
  
}

/* This is the integrand of the streaming model (duh). 
 *
 *                [1+xi(r)]*P(vz)
 *
 */
double streaming_model_integrand(double z)
{
  double r,theta,vx,xi,p,vz;

  r=sqrt(z*z+rs2_g5);
  theta=asin(z/r);
  if(theta==0)theta=1.0E-6; /* Prevent NANs */
  vz=HUBBLEZ*100*(rp_g5-z);

  xi=two_halo_real_space(r);
  if(z<0)vz=-vz; /* This is for interpolation in the histogram table. */
  p=galaxy_prob_vz(vz,r,theta);
  //printf("BOO %f %f %e %e %f\n",z,vz,p,xi,theta/PI*180);
  return(p*(xi+1));
  
}

/* This is a check to make sure that the integral being passed
 * to qromo isn't full of zeros.
 */
int xi2h_zero_check(double z1, double z2)
{
  double z,vz,r,theta,ptot=0,p;
  int i;
  
  for(i=1;i<=5;++i)
    {
      z = (i-1)*(z2-z1)/4+z1;
      r=sqrt(z*z+rs2_g5);
      theta=asin(z/r);
      if(theta==0)theta=1.0E-6; /* Prevent NANs */
      vz=HUBBLEZ*100*(rp_g5-z);
      if(z<0)vz=-vz;
      p=galaxy_prob_vz(vz,r,theta);
      ptot+=p;
    }
  if(ptot>0)return(1);
  return(0);
}
