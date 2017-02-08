#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

struct COSMO {
  double om_mhsq;
  double om_bhsq;
  double h;
  double w;
  double ns;
  double sigma_8;
} cosmo;

int HOD_HYPERCUBE=0, icode_max, PK=1;
double glh_lgm1 = 14.18, glh_dlgm1 = 0.13;
double glh_lgmcut = 13.5, glh_dlgmcut = 0.5;
double glh_alpha = 1.4, glh_dalpha = 0.22;
double glh_slogm = 0.56, glh_dslogm = 0.1;

void random_cosmo_params2();
void hypercube_cosmo_params2();
void random_cosmo_params();
void hypercube_cosmo_params();
void random_hod_params();
void hypercube_hod_params();
void nested_hypercube_hod_params(int icode);
double xi_value1(double a[], double r, double xfid);

void latin_hypercube(int icode)
{
  double klo=0.01, khi=3, dlogk, xk, pkfid[1000], fisher[15][15], deltax[15], a0[15], a[15], x0;
  double chi2xmm[10][10], chi2xpp[10][01], chi2xp[10], chi2xm[10], temp[10][10];
  int i1,i2,i3,i4,i5,nparams = 5, i, nsample=1000, j, nk=100, k, ncosmo;
  FILE *fp, *fp1;
  char fname[1000], aa[1000], command[1000];

  HIGH_PRECISION = 0;

  if(icode<0)
    {
      HOD_HYPERCUBE = 1;
      icode *= -1;
    }

  srand48(555);
  klo = 0.1;
  khi = 40.0;
  nk = 30;
  dlogk = log(khi/klo)/(nk-1);

  if(PK)
    {
      klo = 0.3;
      khi = 10.0;
      nk = 10;
      dlogk = log(khi/klo)/(nk-1);
    }

  // let's get the fiducial cosmology
  SIGMA_8 = 0.75;
  HUBBLE = 0.70;
  OMEGA_M = (0.155 + 0.120)/2./HUBBLE/HUBBLE;
  OMEGA_B = (0.0235 + 0.0215)/2./HUBBLE/HUBBLE;
  SPECTRAL_INDX = 0.95;
  RESET_COSMOLOGY++;

  // let's set the fiducial HOD params (from fit1.fit in HOD_FITS)
  HOD.M_min = 0;
  HOD.M1 = pow(10.0, glh_lgm1);
  HOD.alpha = glh_alpha;
  HOD.M_cut = pow(10.0, glh_lgmcut);
  HOD.sigma_logM = glh_slogm;
  set_HOD_params();
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  if((icode==1 || icode==2) && !HOD_HYPERCUBE)
    {
      fp = fopen("fiducial_model.dat","w");
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = xi_interp(xk);
	  fprintf(fp,"%e %e\n",xk,pkfid[k]);
	}
      fclose(fp);
    }
  else
    {
      fp = fopen("fiducial_model.dat","w");
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = projected_xi(xk);
	  fprintf(fp,"%e %e\n",xk,pkfid[k]);
	}
      fclose(fp);
    }

  if(PK)
    {
      OMEGA_M = 0.3;
      SIGMA_8 = 0.8;
      HUBBLE = 0.7;
      SPECTRAL_INDX = 0.95;
      SPECTRAL_RUN = 0.0;
      GAMMA = OMEGA_M*HUBBLE;
      RESET_COSMOLOGY++;

      fp = fopen("fiducial_model.dat","w");
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = nonlinear_power_spectrum(xk);
	  fprintf(fp,"%e %e\n",xk,pkfid[k]);
	}
      fclose(fp);
    }

  // the grid--- well, how about a bunch of random points?
  if(icode==1) 
    {
      for(i=1; i<=nsample; ++i)
	{
	  random_cosmo_params2();
	  sprintf(fname,"random_point.%d",i);
	  fp = fopen(fname,"w");
	  fprintf(fp,"%e %e %e %e %e ",cosmo.om_mhsq, cosmo.om_bhsq, cosmo.h, cosmo.ns, cosmo.sigma_8);
	  //fprintf(fp,"%e %e %e %e\n",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      fprintf(fp,"%e %e\n",log(xk),log(nonlinear_power_spectrum(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(xi_interp(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	    }
	  fclose(fp);
	}
	  
    }

  // input the hypercube and output the values
  if(icode==2)
    {
      //system("/Users/tinker/src/latin_hypercube/optimal_spacing 5 80 1784391832 1 > lh.out");      // for matter
      //system("/Users/tinker/src/latin_hypercube/optimal_spacing 9 160 555 0 > lh.out");    // for wp  
      //sprintf(fname,"hypercube_input_wide_lhr_1784391832_np5_n80_nr17.dat",i);
      //sprintf(fname,"hypercube_input_wide_lhr_mcmc_2dsum_558_np5_n80_nr17.dat",i);
      sprintf(fname,"hypercube_input_wide_lhr_mcmc_np5_n40.dat",i);
      fp = fopen(fname,"w");
      for(i=1;i<=40;++i)
	{
	  hypercube_cosmo_params2();	  
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      //fmuh(xk);
	      fprintf(fp,"%e %e %e %e %e ",cosmo.om_mhsq, cosmo.om_bhsq, cosmo.h, cosmo.ns, cosmo.sigma_8);
	      //fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	      //fprintf(fp,"%e %e\n",log(xk),log(xi_interp(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	      fprintf(fp,"%e %e\n",log(xk),log(nonlinear_power_spectrum(xk)/pkfid[k]));
	    }
	  fflush(fp);
	}
      fclose(fp);
    }

  // the grid--- well, how about a bunch of random points?
  if(icode==3) 
    {
      for(i=1; i<=nsample; ++i)
	{
	  random_hod_params();
	  sprintf(fname,"random_hod_point.%d",i);
	  fp = fopen(fname,"w");
	  fprintf(fp,"%e %e %e %e\n",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	    }
	  fclose(fp);
	}
	  
    }
  // input the hypercube and output the values
  if(icode==4)
    {
      system("/Users/tinker/src/latin_hypercube/optimal_spacing 4 40 555 0 > lh.out");      
      sprintf(fname,"hypercube_input_wide_lhc_555_np4_n40_hod.dat",i);
      fp = fopen(fname,"w");
      for(i=1;i<=40;++i)
	{
	  hypercube_hod_params();	  
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	      fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	    }
	  fflush(fp);
	}
      fclose(fp);
    }

  // this is for nested hypercubes
  if(icode==5)
    {
      system("/Users/tinker/src/latin_hypercube/optimal_spacing 5 40 557 0 > lh.out");      
      sprintf(fname,"hypercube_input_wide_lhc_557_np5_n40_nested_np4_n20_nr15.dat",i);
      fp = fopen(fname,"w");
      for(i=1;i<=40;++i)
	{
	  hypercube_cosmo_params();
	  icode_max = 20;
	  sprintf(aa,"/Users/tinker/src/latin_hypercube/optimal_spacing 4 %d %d 1 > lhx.out",
		  icode_max,555+i);
	  system(aa);      
	  for(j=1;j<=icode_max;++j)
	    {
	      nested_hypercube_hod_params(j);
	      for(k=15;k<16;++k)
		{
		  xk = klo*exp(k*dlogk);
		  fprintf(fp,"%e %e %e %e %e ",cosmo.om_mhsq, cosmo.om_bhsq, cosmo.h, cosmo.ns, cosmo.sigma_8);
		  fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
		  fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
		}
	      fflush(fp);
	    }
	}
      fclose(fp);
    }
      
  // just a single parameter
  if(icode==6)
    {
      fp = fopen("hypercube_input_short.dat","w");
      for(i=-5;i<=5;++i)
	{
	  GAMMA = 0.28 + 0.04*i/5.;
	  RESET_COSMOLOGY++;
	  fprintf(fp,"%e %e\n",GAMMA,log(nonlinear_power_spectrum(0.3)/pkfid[0]));
	}
      fclose(fp);
      
      for(i=1;i<=1000;++i)
	{
	  GAMMA = 0.28 + 0.04*(drand48()*2-1);
	  RESET_COSMOLOGY++;
	  sprintf(fname,"random_point.%d",i);
	  fp = fopen(fname,"w");
	  fprintf(fp,"%e %e\n",GAMMA,log(nonlinear_power_spectrum(0.3)/pkfid[0]));
	  fclose(fp);
	}
    }


  // inputting hypercube eigen params
  if(icode==7)
    {
      LINEAR_PSPEC_FILE = 1;
      RESET_COSMOLOGY++;

      // lets remake the fiducial cosmology
      klo = 0.3;
      khi = 10.0;
      nk = 10;
      dlogk = log(khi/klo)/(nk-1);

      cosmo.om_bhsq = 2.226292e-02;
      cosmo.om_mhsq = 1.178303e-01; 
      cosmo.w = -9.965125e-01;
      cosmo.ns =  9.625152e-01;
      cosmo.sigma_8 = 8.177711e-01; 
      cosmo.h = 6.823171e+01;

      HUBBLE = cosmo.h = cosmo.h/100;
      OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
      OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
      SIGMA_8 = cosmo.sigma_8;
      SPECTRAL_INDX = cosmo.ns;

      // run the CAMB script
      sprintf(command,"./sh.initfile %e %e %e %e",cosmo.om_bhsq, cosmo.om_mhsq, cosmo.h*100, cosmo.w);
      printf("COMMAND: %s\n",command);
      system(command);
      RESET_COSMOLOGY++;

      fp = fopen("fiducial_model.dat","w");
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = nonlinear_power_spectrum(xk);
	  fprintf(fp,"%e %e\n",xk,pkfid[k]);
	}
      fclose(fp);
      fprintf(stderr,"done with fiducial model\n");

      sprintf(fname,"hypercube_input_wide_n80.dat");
      fp = fopen(fname,"w");
      fp1 = fopen("LH_eigenspace_n80.out","r");      
      ncosmo = filesize(fp1);
      for(i=1;i<=ncosmo;++i)
	{
	  // read in the cosmology
	  fscanf(fp1,"%lf %lf %lf %lf %lf %lf", 
		 &cosmo.om_bhsq, &cosmo.om_mhsq, &cosmo.w, &cosmo.ns, &cosmo.sigma_8, &cosmo.h);
	  HUBBLE = cosmo.h = cosmo.h/100;
	  OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
	  OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
	  SIGMA_8 = cosmo.sigma_8;
	  SPECTRAL_INDX = cosmo.ns;

	  // run the CAMB script
	  sprintf(command,"./sh.initfile %e %e %e %e",cosmo.om_bhsq, cosmo.om_mhsq, cosmo.h*100, cosmo.w);
	  printf("COMMAND: %s\n",command);
	  system(command);
	  RESET_COSMOLOGY++;
    
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      //fmuh(xk);
	      fprintf(fp,"%e %e %e %e %e %e ",cosmo.om_bhsq, cosmo.om_mhsq, 
		      cosmo.w, cosmo.ns, cosmo.sigma_8, cosmo.h);
	      //fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	      //fprintf(fp,"%e %e\n",log(xk),log(xi_interp(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	      fprintf(fp,"%e %e\n",log(xk),log(nonlinear_power_spectrum(xk)/pkfid[k]));
	    }
	  fflush(fp);
	}
      fclose(fp);
    }

  // inputting random params
  if(icode==8)
    {
      LINEAR_PSPEC_FILE = 1;
      RESET_COSMOLOGY++;

      // lets remake the fiducial cosmology
      klo = 0.3;
      khi = 10.0;
      nk = 10;
      dlogk = log(khi/klo)/(nk-1);

      cosmo.om_bhsq = 2.226292e-02;
      cosmo.om_mhsq = 1.178303e-01; 
      cosmo.w = -9.965125e-01;
      cosmo.ns =  9.625152e-01;
      cosmo.sigma_8 = 8.177711e-01; 
      cosmo.h = 6.823171e+01;

      HUBBLE = cosmo.h = cosmo.h/100;
      OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
      OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
      SIGMA_8 = cosmo.sigma_8;
      SPECTRAL_INDX = cosmo.ns;

      // run the CAMB script
      sprintf(command,"./sh.initfile %e %e %e %e",cosmo.om_bhsq, cosmo.om_mhsq, cosmo.h*100, cosmo.w);
      printf("COMMAND: %s\n",command);
      system(command);
      RESET_COSMOLOGY++;

      fp = fopen("fiducial_model.dat","w");
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = nonlinear_power_spectrum(xk);
	  fprintf(fp,"%e %e\n",xk,pkfid[k]);
	}
      fclose(fp);
      fprintf(stderr,"done with fiducial model\n");

      fp1 = fopen("sampled_eigenspace.dat","r");      
      for(i=1;i<=1000;++i)
	{
	  sprintf(fname,"random_point.%d",i);
	  fp = fopen(fname,"w");

	  // read in the cosmology
	  fscanf(fp1,"%lf %lf %lf %lf %lf %lf", 
		 &cosmo.om_bhsq, &cosmo.om_mhsq, &cosmo.w, &cosmo.ns, &cosmo.sigma_8, &cosmo.h);
	  HUBBLE = cosmo.h = cosmo.h/100;
	  OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
	  OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
	  SIGMA_8 = cosmo.sigma_8;
	  SPECTRAL_INDX = cosmo.ns;

	  // run the CAMB script
	  sprintf(command,"./sh.initfile %e %e %e %e",cosmo.om_bhsq, cosmo.om_mhsq, cosmo.h*100, cosmo.w);
	  printf("COMMAND: %s\n",command);
	  system(command);
	  RESET_COSMOLOGY++;

	  fprintf(fp,"%e %e %e %e %e %e ",cosmo.om_bhsq, cosmo.om_mhsq, 
		  cosmo.w, cosmo.ns, cosmo.sigma_8, cosmo.h);
    
	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      //fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	      //fprintf(fp,"%e %e\n",log(xk),log(xi_interp(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	      fprintf(fp,"%e %e\n",log(xk),log(nonlinear_power_spectrum(xk)/pkfid[k]));
	    }
	  fflush(fp);
	}
      fclose(fp);
    }

  if(icode==9)
    {
      LINEAR_PSPEC_FILE = 1;
      RESET_COSMOLOGY++;

      // lets remake the fiducial cosmology
      klo = 0.3;
      khi = 10.0;
      nk = 10;
      dlogk = log(khi/klo)/(nk-1);

      cosmo.om_bhsq = 2.226292e-02;
      cosmo.om_mhsq = 1.178303e-01; 
      cosmo.w = -9.965125e-01;
      cosmo.ns =  9.625152e-01;
      cosmo.sigma_8 = 8.177711e-01; 
      cosmo.h = 6.823171e+01;

      HUBBLE = cosmo.h = cosmo.h/100;
      OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
      OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
      SIGMA_8 = cosmo.sigma_8;
      SPECTRAL_INDX = cosmo.ns;
      // run the CAMB script
      sprintf(command,"./sh.initfile %e %e %e %e %d",cosmo.om_bhsq, cosmo.om_mhsq, cosmo.h*100, cosmo.w, 0);
      printf("COMMAND: %s\n",command);
      system(command);
      RESET_COSMOLOGY++;

      printf("NORM fid %e\n",OMEGA_M);
      for(k=0;k<nk;++k)
	{
	  xk = klo*exp(k*dlogk);
	  pkfid[k] = nonlinear_power_spectrum(xk);
	}
      fprintf(stderr,"done with fiducial model\n");

      // lets step through OMEGA_M
      for(i=1;i<=100;++i)
	{
	  sprintf(fname,"omega2_point.%d",i);
	  fp = fopen(fname,"w");

	  OMEGA_M = (0.2*i/100+0.9)*cosmo.om_mhsq/HUBBLE/HUBBLE;

	  // run the CAMB script
	  sprintf(command,"./sh.initfile %e %e %e %e %d",cosmo.om_bhsq, 
		  cosmo.om_mhsq*(0.2*i/100+0.9), cosmo.h*100, cosmo.w, i);
	  printf("COMMAND: %s\n",command);
	  system(command);
	  RESET_COSMOLOGY++;

	  fprintf(fp,"%e %e %e %e %e %e\n ",cosmo.om_bhsq, cosmo.om_mhsq, 
		  cosmo.w, cosmo.ns, cosmo.sigma_8, cosmo.h);
	  fprintf(stderr,"SIGMA_8: %f\n",cosmo.sigma_8);

	  for(k=0;k<nk;++k)
	    {
	      xk = klo*exp(k*dlogk);
	      //fprintf(fp,"%e %e %e %e ",log10(HOD.M1), HOD.alpha, log10(HOD.M_cut), HOD.sigma_logM);
	      //fprintf(fp,"%e %e\n",log(xk),log(xi_interp(xk)/pkfid[k]));
	      //fprintf(fp,"%e %e\n",log(xk),log(projected_xi(xk)/pkfid[k]));
	      fprintf(fp,"%e %e\n",log(xk),log(nonlinear_power_spectrum(xk)/pkfid[k]));
	    }
	  fflush(fp);
	}
      

    }

}

double xi_value1(double a[], double r, double xfid)
{
  SIGMA_8 = a[5];
  HUBBLE = a[3];
  OMEGA_M = a[1]/HUBBLE/HUBBLE;
  OMEGA_B = a[2]/HUBBLE/HUBBLE;
  SPECTRAL_INDX = a[4];
  RESET_COSMOLOGY++;

  return log(xi_interp(r)/xfid);
}

void hypercube_cosmo_params2()
{
  static int iflag = 1;
  static FILE *fp;
  float xx;
  double dom, ds8, dh0, dns, drun;

  if(iflag) {
    fp = openfile("lh.out");
    iflag = 0;
  }
  
  // central model
  OMEGA_M = 0.3;
  SIGMA_8 = 0.8;
  HUBBLE = 0.7;
  SPECTRAL_INDX = 0.95;
  SPECTRAL_RUN = 0.0;

  // wide withs;
  dom = 0.05*1.1;
  ds8 = 0.15*1.1;
  dh0 = 0.05*1.1;
  dns = 0.05*1.1;
  drun = 0.02*1.1;

  fscanf(fp,"%f",&xx);
  OMEGA_M = OMEGA_M + (2*xx-1)*dom;
  fscanf(fp,"%f",&xx);
  SIGMA_8 = SIGMA_8 + (2*xx-1)*ds8;
  fscanf(fp,"%f",&xx);
  HUBBLE = HUBBLE + (2*xx-1)*dh0;
  fscanf(fp,"%f",&xx);
  SPECTRAL_INDX = SPECTRAL_INDX + (2*xx-1)*dns;
  fscanf(fp,"%f",&xx);
  SPECTRAL_RUN = SPECTRAL_RUN + (2*xx-1)*drun;

  GAMMA = HUBBLE*OMEGA_M;
  RESET_COSMOLOGY++;

  cosmo.om_mhsq = OMEGA_M;
  cosmo.om_bhsq = SIGMA_8;
  cosmo.h = HUBBLE;
  cosmo.ns = SPECTRAL_INDX; 
  cosmo.sigma_8 = SPECTRAL_RUN;

}

void random_cosmo_params2()
{
  static int iflag = 1;
  static FILE *fp;
  float xx;
  double dom, ds8, dh0, dns, drun;


  // central model
  OMEGA_M = 0.3;
  SIGMA_8 = 0.8;
  HUBBLE = 0.7;
  SPECTRAL_INDX = 0.95;
  SPECTRAL_RUN = 0.0;

  // wide withs;
  dom = 0.05*1.;
  ds8 = 0.15*1.;
  dh0 = 0.05*1.;
  dns = 0.05*1.;
  drun = 0.02*1.;

  xx = drand48();
  OMEGA_M = OMEGA_M + (2*xx-1)*dom;
  xx = drand48();
  SIGMA_8 = SIGMA_8 + (2*xx-1)*ds8;
  xx = drand48();
  HUBBLE = HUBBLE + (2*xx-1)*dh0;
  xx = drand48();
  SPECTRAL_INDX = SPECTRAL_INDX + (2*xx-1)*dns;
  xx = drand48();
  SPECTRAL_RUN = SPECTRAL_RUN + (2*xx-1)*drun;

  GAMMA = HUBBLE*OMEGA_M;
  RESET_COSMOLOGY++;

  cosmo.om_mhsq = OMEGA_M;
  cosmo.om_bhsq = SIGMA_8;
  cosmo.h = HUBBLE;
  cosmo.ns = SPECTRAL_INDX; 
  cosmo.sigma_8 = SPECTRAL_RUN;
}

void hypercube_cosmo_params()
{
  static int iflag = 1;
  static FILE *fp;
  float xx;
  double sigma8_lo, sigma8_hi, omh2_lo, omh2_hi, ombh2_lo, ombh2_hi, ns_lo, ns_hi, w_hi, w_lo, omnu_hi, h_lo, h_hi; 

  double m1, mcut, alpha, slogm, dm1, dalpha, dmcut, dslogm;

  if(iflag) {
    fp = openfile("lh.out");
    iflag = 0;
  }

  sigma8_lo = 0.6;
  sigma8_hi = 0.9;
  omh2_lo = 0.120;
  omh2_hi = 0.155;
  ombh2_lo = 0.0215;
  ombh2_hi = 0.0235;
  ns_lo = 0.85;
  ns_hi = 1.05;
  h_lo = 0.65;
  h_hi = 0.75;

  // "wide" input
  sigma8_lo = 0.55;
  sigma8_hi = 0.95;
  omh2_lo = 0.115;
  omh2_hi = 0.160;
  ombh2_lo = 0.0212;
  ombh2_hi = 0.0238;
  ns_lo = 0.82;
  ns_hi = 1.08;
  h_lo = 0.62;
  h_hi = 0.78;


  fscanf(fp,"%f",&xx);
  cosmo.om_mhsq = omh2_lo + xx*(omh2_hi - omh2_lo);
  fscanf(fp,"%f",&xx);
  cosmo.om_bhsq = ombh2_lo + xx*(ombh2_hi - ombh2_lo);
  fscanf(fp,"%f",&xx);
  cosmo.ns = ns_lo + xx*(ns_hi - ns_lo);
  fscanf(fp,"%f",&xx);
  cosmo.h = h_lo + xx*(h_hi - h_lo);
  fscanf(fp,"%f",&xx);
  cosmo.sigma_8 = sigma8_lo + xx*(sigma8_hi - sigma8_lo);

  HUBBLE = cosmo.h;
  OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
  OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
  SPECTRAL_INDX = cosmo.ns;
  SIGMA_8 = cosmo.sigma_8;
  RESET_COSMOLOGY++;
  
  // addint the 1.25 to each "d" for the "wide" input.
  if(HOD_HYPERCUBE)
    {
      m1 = log10(1.272076e+14);
      alpha = 1.49;
      mcut = log10(5.0E13);
      slogm = 0.54;
      dm1 = 0.04*1.25;
      dalpha = 0.15*1.25;
      dmcut = 0.5*1.25;
      dslogm = 0.05*1.25;

      fscanf(fp,"%f",&xx);
      HOD.M1 = pow(10.0,m1 + 2*(xx-0.5)*dm1);
      fscanf(fp,"%f",&xx);
      HOD.alpha = alpha + 2*(xx-0.5)*dalpha;
      fscanf(fp,"%f",&xx);
      HOD.M_cut = pow(10.0,mcut + 2*(xx-0.5)*dmcut);
      fscanf(fp,"%f",&xx);
      HOD.sigma_logM = slogm + 2*(xx-0.5)*dslogm;

      HOD.M_min = 0;
      set_HOD_params();
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;

    }

}
void hypercube_hod_params()
{
  static int iflag = 1;
  static FILE *fp;
  float xx, dfac=1.15;
  double sigma8_lo, sigma8_hi, omh2_lo, omh2_hi, ombh2_lo, ombh2_hi, ns_lo, ns_hi, w_hi, w_lo, omnu_hi, h_lo, h_hi; 
  double m1, mcut, alpha, slogm, dm1, dalpha, dmcut, dslogm;

  if(iflag) {
    fp = openfile("lh.out");
    iflag = 0;
  }

  fscanf(fp,"%f",&xx);
  HOD.M1 = pow(10.0,glh_lgm1 + 2*(xx-0.5)*glh_dlgm1*dfac);
  fscanf(fp,"%f",&xx);
  HOD.alpha = glh_alpha + 2*(xx-0.5)*glh_dalpha*dfac;
  fscanf(fp,"%f",&xx);
  HOD.M_cut = pow(10.0,glh_lgmcut + 2*(xx-0.5)*glh_dlgmcut*dfac);
  fscanf(fp,"%f",&xx);
  HOD.sigma_logM = glh_slogm + 2*(xx-0.5)*glh_dslogm*dfac;
  
  HOD.M_min = 0;
  set_HOD_params();
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;
  
}
void nested_hypercube_hod_params(int icode)
{
  static int iflag = 1;
  static FILE *fp;
  float xx, dfac=1.15;
  double sigma8_lo, sigma8_hi, omh2_lo, omh2_hi, ombh2_lo, ombh2_hi, ns_lo, ns_hi, w_hi, w_lo, omnu_hi, h_lo, h_hi; 
  double m1, mcut, alpha, slogm, dm1, dalpha, dmcut, dslogm;

  if(icode==1) {
    fp = openfile("lhx.out");
  }

  fscanf(fp,"%f",&xx);
  HOD.M1 = pow(10.0,glh_lgm1 + 2*(xx-0.5)*glh_dlgm1*dfac);
  fscanf(fp,"%f",&xx);
  HOD.alpha = glh_alpha + 2*(xx-0.5)*glh_dalpha*dfac;
  fscanf(fp,"%f",&xx);
  HOD.M_cut = pow(10.0,glh_lgmcut + 2*(xx-0.5)*glh_dlgmcut*dfac);
  fscanf(fp,"%f",&xx);
  HOD.sigma_logM = glh_slogm + 2*(xx-0.5)*glh_dslogm*dfac;
  
  HOD.M_min = 0;
  set_HOD_params();
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;
  
  if(icode==icode_max) {
    fclose(fp);
  }

}

void random_cosmo_params()
{
  double sigma8_lo, sigma8_hi, omh2_lo, omh2_hi, ombh2_lo, ombh2_hi, ns_lo, ns_hi, w_hi, w_lo, omnu_hi, h_lo, h_hi; 
  double m1, mcut, alpha, slogm, dm1, dalpha, dmcut, dslogm, xx;

  // same as the hypercube (not the "wide")
  sigma8_lo = 0.6;
  sigma8_hi = 0.9;
  omh2_lo = 0.120;
  omh2_hi = 0.155;
  ombh2_lo = 0.0215;
  ombh2_hi = 0.0235;
  ns_lo = 0.85;
  ns_hi = 1.05;
  h_lo = 0.65;
  h_hi = 0.75;

  // "narrow"
  /*
  sigma8_lo = 0.65;
  sigma8_hi = 0.85;
  omh2_lo = 0.125;
  omh2_hi = 0.150;
  ombh2_lo = 0.0218;
  ombh2_hi = 0.0232;
  ns_lo = 0.88;
  ns_hi = 1.02;
  h_lo = 0.68;
  h_hi = 0.72;
  */


  cosmo.om_mhsq = omh2_lo + drand48()*(omh2_hi - omh2_lo);
  cosmo.om_bhsq = ombh2_lo + drand48()*(ombh2_hi - ombh2_lo);
  cosmo.ns = ns_lo + drand48()*(ns_hi - ns_lo);
  cosmo.h = h_lo + drand48()*(h_hi - h_lo);
  cosmo.sigma_8 = sigma8_lo + drand48()*(sigma8_hi - sigma8_lo);

  HUBBLE = cosmo.h;
  OMEGA_M = cosmo.om_mhsq/HUBBLE/HUBBLE;
  OMEGA_B = cosmo.om_bhsq/HUBBLE/HUBBLE;
  SPECTRAL_INDX = cosmo.ns;
  SIGMA_8 = cosmo.sigma_8;
  RESET_COSMOLOGY++;

  if(HOD_HYPERCUBE)
    {
      xx = drand48();
      HOD.M1 = pow(10.0,glh_lgm1 + 2*(xx-0.5)*glh_dlgm1);
      xx = drand48();
      HOD.alpha = glh_alpha + 2*(xx-0.5)*glh_dalpha;
      xx = drand48();
      HOD.M_cut = pow(10.0,glh_lgmcut + 2*(xx-0.5)*glh_dlgmcut);
      xx = drand48();
      HOD.sigma_logM = glh_slogm + 2*(xx-0.5)*glh_dslogm;

      HOD.M_min = 0;
      set_HOD_params();
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;

    }

}
void random_hod_params()
{
  double sigma8_lo, sigma8_hi, omh2_lo, omh2_hi, ombh2_lo, ombh2_hi, ns_lo, ns_hi, w_hi, w_lo, omnu_hi, h_lo, h_hi; 
  double m1, mcut, alpha, slogm, dm1, dalpha, dmcut, dslogm, xx;

  xx = drand48();
  HOD.M1 = pow(10.0,glh_lgm1 + 2*(xx-0.5)*glh_dlgm1);
  xx = drand48();
  HOD.alpha = glh_alpha + 2*(xx-0.5)*glh_dalpha;
  xx = drand48();
  HOD.M_cut = pow(10.0,glh_lgmcut + 2*(xx-0.5)*glh_dlgmcut);
  xx = drand48();
  HOD.sigma_logM = glh_slogm + 2*(xx-0.5)*glh_dslogm;

  HOD.M_min = 0;
  set_HOD_params();
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

}

double  func_cosmo_params(double h)
{
  double om_mhsq, om_bhsq, e, thet, thetsq, thetpf, b1, b2, zd, za, rd, re, rke, s, ze, dls;

  //om_mhsq = cosmo.om_mhsq;
  //om_bhsq = cosmo.om_bhsq;

  // constants
  e=exp(1.);      
  thet=2.728/2.7;
  thetsq=thet*thet;
  thetpf=thetsq*thetsq;

  // Equation 4 - redshift of drag epoch
  b1=0.313*pow(om_mhsq,-0.419)*(1.+0.607*pow(om_mhsq,0.674));
  b2=0.238*pow(om_mhsq,0.223);
  zd=1291.*(1.+b1*pow(om_bhsq,b2))*pow(om_mhsq,0.251)
    /(1.+0.659*pow(om_mhsq,0.828));

  // Equation 2 - redshift of matter-radiation equality
  ze=2.50e4*om_mhsq/thetpf;

  // value of R=(ratio of baryon-photon momentum density) at drag epoch
  rd=31500.*om_bhsq/(thetpf*zd);

  // value of R=(ratio of baryon-photon momentum density) at epoch of matter-radiation equality
  re=31500.*om_bhsq/(thetpf*ze);

  // Equation 3 - scale of ptcle horizon at matter-radiation equality
  rke=7.46e-2*om_mhsq/(thetsq);

  // Equation 6 - sound horizon at drag epoch
  s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re)));

  // now solve for distance to the drag epoch (requires h)
  OMEGA_M = om_mhsq/h/h;  
  dls = distance_redshift(zd);
  return PI*dls/rd - 302.4; // use the constraint that ell_A = 302.4

}
