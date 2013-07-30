#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"

void cisi(double x, double *ci, double *si);

/* The Fourier transform of the NFW profile.
 * Inputs are wavenumber and halo mass.
 * NB! -> Function assumes that profile of satellite GALAXIES is
 * desired, therefore uses CVIR_FAC
 */

double nfw_transform(double xk, double m)
{
  double c,rvir,kappa1,kappa2,f,y,ci1,ci2,si1,si2;

  c=halo_concentration(m)*CVIR_FAC;
  rvir=pow(3.0*m/(4.0*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  kappa1=xk*rvir/c;
  kappa2=kappa1*(1.0+c);
  f=log(1.0+c)-c/(1.0+c);
  f=1.0/f;
  
  cisi(kappa1,&ci1,&si1);
  cisi(kappa2,&ci2,&si2);

  y=-sin(kappa1*c)/kappa2;
  y=y+sin(kappa1)*(si2-si1)+cos(kappa1)*(ci2-ci1);
  y=f*y;

  return(y);
}

