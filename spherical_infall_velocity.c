#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"

double find_theta(double theta);
double delta1;

/* NOTE--> The halo mass must be in units of h^{-1} M_sol 
 */
/* NOTE--> Wed Jan 28 10:18:13 EST 2004 Tyring reducing the halo mass used in the
 * calculation (was 2*halo_mass, now 1*halo_mass) > REMOVED
 */

/* This is a different version. Just give it delta and it will give you the velocity.
 * (Independent of radius.)
 */
double spherical_collapse_model(double delta)
{
  double theta,u;

  delta1=delta;
  theta=zbrent(find_theta,0.0001,1.99*PI,1.0E-4);  
  u=-(3.0*sin(theta)*(theta-sin(theta))/(2*(1-cos(theta))*(1-cos(theta)))-1);
  u=-u*100*pow(OMEGA_M,0.6);
  return(u);
}

double find_theta(double theta)
{
  return(4.5*(theta-sin(theta))*(theta-sin(theta))/(pow(1-cos(theta),3.0))-1-delta1);
}

