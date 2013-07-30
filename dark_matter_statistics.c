#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"

void output_matter_power_spectrum()
{
  int k;
  double dlogk,xk,plin,pnl;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING MATTER POWER SPECTRA.\n");
  fprintf(stderr,    "---------------------------------\n\n");

  sprintf(aa,"%s.matter_pk",Task.root_filename);
  fp = fopen(aa,"w");

  for(k=-300;k<=100;++k)
    {
      xk=pow(10.0,k/100.0);
      plin = linear_power_spectrum(xk)/(xk*xk*xk)*TWOPI*PI;
      pnl = nonlinear_power_spectrum(xk)/(xk*xk*xk)*TWOPI*PI;
      fprintf(fp,"%e %e %e\n",xk,plin,pnl);
    }
  fclose(fp);

}

void output_matter_correlation_function()
{
  int k,nr=50;
  double dlogr,r,xlin,xnl,rmin = 0.05,rmax = 80.0;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING MATTER CORRELATION FUNCTIONIONS.\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  sprintf(aa,"%s.matter_xi",Task.root_filename);
  fp = fopen(aa,"w");

  dlogr = (log(rmax) - log(rmin))/(nr-1);

  for(k=0;k<nr;++k)
    {
      r = exp(k*dlogr)*rmin;
      xlin = xi_linear_interp(r);
      xnl = xi_interp(r);
      fprintf(fp,"%e %e %e %e\n",r,xlin,xnl,projected_xi_matter(r));
    }
  fclose(fp);
}

void output_halo_correlation_function(double mass)
{
  int k,nr=50;
  double dlogr,r,xlin,xnl,rmin = 0.15,rmax = 40.0;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO CORRELATION FUNCTIONION (M=%.2e)\n",mass);
  fprintf(stderr,    "-------------------------------------------------\n\n");

  sprintf(aa,"%s.halo_xi",Task.root_filename);
  fp = fopen(aa,"w");

  rmin = pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),THIRD);

  dlogr = (log(rmax) - log(rmin))/(nr-1);
  
  for(k=0;k<nr;++k)
    {
      r = exp(k*dlogr)*rmin;
      xlin =  bias_interp(mass,r); 
      xnl = xi_interp(r);
      fprintf(fp,"%e %e\n",r,xnl*xlin*xlin);
    }
  fclose(fp);
}

void output_matter_variance()
{
  int k,nr=50;
  double dlogr,r,slin,snl,mass,rmin = 0.05,rmax = 80.0,pnorm,pnorm_nl;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING MATTER VARIANCE.\n");
  fprintf(stderr,    "----------------------------\n\n");

  sprintf(aa,"%s.sigma_r",Task.root_filename);
  fp = fopen(aa,"w");

  dlogr = (log(rmax) - log(rmin))/(nr-1);
  pnorm = SIGMA_8/sigmac(8.0);
  pnorm_nl = pnorm*sigmac(80.0)/nonlinear_sigmac(80.0);

  for(k=0;k<nr;++k)
    {
      r = exp(k*dlogr)*rmin;
      mass = 4./3.*PI*r*r*r*RHO_CRIT*OMEGA_M;
      slin = sigmac(r)*pnorm;
      snl = nonlinear_sigmac(r)*pnorm_nl;
      fprintf(fp,"%e %e %e %e\n",r,slin,snl,mass);
    }
  fclose(fp);

}

void output_halo_concentrations()
{
  int k,nr=100;
  double x,dlogm,m,mvir,cdelta,mmin=1e8,mmax=1e16,delta_vir;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO CONCENTRATIONS.\n");
  fprintf(stderr,    "--------------------------------\n\n");

  sprintf(aa,"%s.cvir",Task.root_filename);
  fp = fopen(aa,"w");

  dlogm = (log(mmax) - log(mmin))/(nr-1);

  for(k=0;k<nr;++k)
    {
      m = exp(k*dlogm)*mmin;
      cdelta = halo_concentration(m);
      fprintf(fp,"%e %e\n",m,cdelta);
    }
  fclose(fp);

}

void output_halo_mass_function()
{
  int k,nr=100;
  double x,dlogm,m,mvir,cdelta,mmin=1e8,mmax=1e16,delta_vir,b;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO MASS FUNCTION + BIAS FUNCTION.\n");
  fprintf(stderr,    "----------------------------------------------\n\n");

  sprintf(aa,"%s.halostats",Task.root_filename);
  fp = fopen(aa,"w");

  dlogm = (log(mmax) - log(mmin))/(nr-1);

  for(k=0;k<nr;++k)
    {
      m = exp(k*dlogm)*mmin;
      x = dndM_interp(m);
      b = bias(m);
      fprintf(fp,"%e %e %e\n",m,x,b);
    }
  fclose(fp);
}
