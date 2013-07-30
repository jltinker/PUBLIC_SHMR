#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"


/* This is a routine which takes a virial mass (as defined by the virial
 * overdensity relation (see fitting formula of Bryan & Normax)
 * and converts if to a mass of some given overdensity. This is all
 * taken from the appendix of Hu & Kravtsov. 2003ApJ...584..702H
 */

/* Constants from conversion formula in Hu & Kravtsov (2003)
 */
#define a1 0.5116
#define a2 -0.4283
#define a3 -3.13E-3
#define a4 -3.52E-5

double HK_func(double x);

double halo_mass_conversion(double mvir, double *cvir1, double delta_halo)
{
  double x,y,z,vx,vy,vz,x1,x2,x3,x4,omega_m,delta=180,mstar,delta_vir,delta_fof,p,f,
    mp,v,mass,lambda,error,tol=1.0E-3,cprev,m_delta,delta_fac,cvir,rdelta,rvir;

  cvir=*cvir1;
  omega_m=OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  lambda=1-OMEGA_M;
  mstar=MSTAR;
  delta=delta_halo;

  x=omega_m-1;
  delta_vir=(18*PI*PI+82*x-39*x*x)/(1+x);

  /* Now convert from virial mass to m_delta.
   */
  f=delta/delta_vir*HK_func(1.0/cvir);
  p=a2+a3*log(f)+a4*log(f)*log(f);
  x1=1.0/sqrt(a1*pow(f,2*p)+0.5625)+2*f;
  m_delta=mvir*(delta/delta_vir*pow(1.0/(x1*cvir),3.0));

  /* Now convert the cvir to c_delta.
   */
  rvir=pow(3*mvir/(4*delta_vir*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rdelta=pow(3*m_delta/(4*delta*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  *cvir1=cvir*rdelta/rvir;

  return(m_delta);
}

/* This is a slight modification to the above routine-- instead of converting from
 * the virial mass, it converts from a specified Delta to the other Delta. 
 * (so here the "_vir" quantities are any arbitrary input overdensity.)
 */

double halo_mass_conversion2(double mvir, double cvir1, double delta_vir, double delta_halo)
{
  double x,y,z,vx,vy,vz,x1,x2,x3,x4,omega_m,delta=180,mstar,delta_fof,p,f,
    mp,v,mass,lambda,error,tol=1.0E-3,cprev,m_delta,delta_fac,cvir,rdelta,rvir;

  cvir=cvir1;
  omega_m=OMEGA_M;
  lambda=1-OMEGA_M;
  mstar=MSTAR;
  delta=delta_halo;

  /* Now convert from virial mass to m_delta.
   */
  f=delta/delta_vir*HK_func(1.0/cvir);
  p=a2+a3*log(f)+a4*log(f)*log(f);
  x1=1.0/sqrt(a1*pow(f,2*p)+0.5625)+2*f;
  m_delta=mvir*(delta/delta_vir*pow(1.0/(x1*cvir),3.0));

  /* Now convert the cvir to c_delta.
   */
  rvir=pow(3*mvir/(4*delta_vir*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rdelta=pow(3*m_delta/(4*delta*PI*RHO_CRIT*OMEGA_M),1.0/3.0);

  return(m_delta);
}

double HK_func(double x)
{
  return(x*x*x*(log(1+1.0/x)-1.0/(1+x)));
}
