#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* ------------------------------------------------
 * THIS IS TINKER's VERSION OF CHI2_LENSING:
 * DON'T FORGET: INPUT IS IN PHYSICAL H=72 UNITS!!
 * ------------------------------------------------
 */

// local functions
void lensing_data_input(void);
void set_lensing_mass_bins(void);

// local global variables
int LENSING_COVAR = 1;

double chi2_lensing(double *a)
{
  static int flag = 1, iter = 0;
  double mlo, mhi, r, rlo, rhi, dlogr, chi2=0, chi2red, chi2tot;
  int i, j, k, nr;
  FILE *fp;
  char fname[1000];

  nr = 100;
  rlo = 0.01;
  rhi = 10.0;
  dlogr = log(rhi/rlo)/(nr-1);

  if(flag)
    lensing_data_input();

  if(!flag)
    wpx.calculate_two_halo = 0;

  if(LENSING_OUTPUT_FLAG)
    {
      sprintf(fname,"lensing_outfile_red.z%d",wpl.iz);
      fp = fopen(fname,"w");
    }

  chi2tot = 0;
  BLUE_FLAG = 0;
  wpl.reset_inversion = 1;
  for(j=1;j<=wpx.nsamples;++j)
    {
      mlo = pow(10.0,wpx.mstar_wplo[2][j]);
      mhi = pow(10.0,wpx.mstar_wphi[2][j]);
      set_up_hod_for_shmr(mlo,mhi,a);
      
      // just output the xigm
      /*
      delta_sigma(1);
      for(i=1;i<=-nr;++i)
	{
	  r = exp((i-1)*dlogr)*rlo;
	  printf("XX 2 %d %e %e %e\n",j,r,delta_sigma(r),one_halo_galaxy_matter(r)+two_halo_real_space(r));
	}
      */
      for(i=1;i<=wpx.ndatar[j];++i)
	{
	  r = wpx.rdatar[j][i];
	  wpx.model[j][i] = delta_sigma(r);
	  if(!wpx.calculate_two_halo) wpx.model[j][i] += wpx.two_halo_red[j][i];
	  chi2 += (wpx.model[j][i] - wpx.xdatar[j][i])*(wpx.model[j][i] - wpx.xdatar[j][i])/
	    (wpx.edatar[j][i]*wpx.edatar[j][i]);
	  if(LENSING_OUTPUT_FLAG)
	    {
	      fprintf(fp,"%d %e %e %e %e %e\n",j,r*1000.0/HUBBLE/(1+REDSHIFT),
		      wpx.xdatar[j][i],
		      wpx.edatar[j][i],delta_sigma(r),one_halo_galaxy_matter(r));
	    }
	}
      if(flag)
	{
	  wpx.calculate_two_halo = 0;
	  for(i=1;i<=wpx.ndatar[j];++i)
	    {
	      r = wpx.rdatar[j][i];
	      wpx.two_halo_red[j][i] = wpx.model[j][i] - delta_sigma(r);
	    }
	  wpx.calculate_two_halo = 1;
	}	
      RESET_FLAG_DS = 1;
      if(LENSING_COVAR)
	{
	  chi2 = 0;
	  for(i=1;i<=wpx.ndatar[j];++i)
	    for(k=1;k<=wpx.ndatar[j];++k)
	      chi2 += (wpx.model[j][i]-wpx.xdatar[j][i])*(wpx.model[j][k]-wpx.xdatar[j][k])*
		wpx.covarr[j][i][k];
	  chi2tot += chi2;
	}
    }
  if(LENSING_OUTPUT_FLAG)
    fclose(fp);
  chi2red = chi2tot;

  if(LENSING_OUTPUT_FLAG)
    {
      sprintf(fname,"lensing_outfile_blue.z%d",wpl.iz);
      fp = fopen(fname,"w");
    }

  chi2tot = 0;
  BLUE_FLAG = 1;
  wpl.reset_inversion = 1;
  for(j=1;j<=wpx.nsamples;++j)
    {
      mlo = pow(10.0,wpx.mstar_wplo[1][j]);
      mhi = pow(10.0,wpx.mstar_wphi[1][j]);
      set_up_hod_for_shmr(mlo,mhi,a);
      
      for(i=1;i<=wpx.ndatab[j];++i)
	{
	  r = wpx.rdatab[j][i];
	  wpx.model[j][i] = delta_sigma(r);
	  if(!wpx.calculate_two_halo) wpx.model[j][i] += wpx.two_halo_blue[j][i];
	  chi2 += (wpx.model[j][i] - wpx.xdatab[j][i])*(wpx.model[j][i] - wpx.xdatab[j][i])/
	    (wpx.edatab[j][i]*wpx.edatab[j][i]);
	  if(LENSING_OUTPUT_FLAG)
	    {
	      fprintf(fp,"%d %e %e %e %e %e\n",j,r*1000.0/HUBBLE/(1+REDSHIFT),
		      wpx.xdatab[j][i],
		      wpx.edatab[j][i],delta_sigma(r),one_halo_galaxy_matter(r));
	    }
	}
      if(flag)
	{
	  wpx.calculate_two_halo = 0;
	  for(i=1;i<=wpx.ndatab[j];++i)
	    {
	      r = wpx.rdatar[j][i];
	      wpx.two_halo_blue[j][i] = wpx.model[j][i] - delta_sigma(r);
	    }
	  wpx.calculate_two_halo = 1;
	}	
      RESET_FLAG_DS = 1;
      if(LENSING_COVAR)
	{
	  chi2 = 0;
	  for(i=1;i<=wpx.ndatab[j];++i)
	    for(k=1;k<=wpx.ndatab[j];++k)
	      chi2 += (wpx.model[j][i]-wpx.xdatab[j][i])*(wpx.model[j][k]-wpx.xdatab[j][k])*
		wpx.covarb[j][i][k];
	  chi2tot += chi2;
	}
    }
  if(LENSING_OUTPUT_FLAG)
    fclose(fp);
  if(LENSING_COVAR)chi2 = chi2tot + chi2red;

  if(flag)flag = 0;
  printf("CHILENSING %d %e %e %e\n",iter++,chi2,chi2red,chi2-chi2red);
  return chi2;
}


void lensing_data_input(void)
{
  double x1;
  int i,j,n,k,i1, istart, ii, nlim = 20;
  FILE *fp;
  char fname[1000], aa[1000], color[100];
  double **covar, **tmp, **tmp2;

  wpx.zlo = wpx.zhi = -1;
  if(REDSHIFT>=0.2 || REDSHIFT<0.5)
    {
      wpx.iz = 1;
      wpx.zlo = 0.22;
      wpx.zhi = 0.48;
      wpx.nsamples = 7;
    }
  if(REDSHIFT>=0.5 && REDSHIFT <0.75)
    {
      wpx.iz = 2;
      wpx.zlo = 0.48;
      wpx.zhi = 0.74;
      wpx.nsamples = 7;
    }
  if(REDSHIFT>=0.75 && REDSHIFT <1.00)
    {
      wpx.iz = 3;
      wpx.zlo = 0.74;
      wpx.zhi = 1.00;
      wpx.nsamples = 6;
    }
  istart = 1;
  if(wpx.zlo<0)
    endrun("ERROR: no proper redshift specified.\n");


  wpx.mstar_threshold = dvector(1,wpx.nsamples);
  wpx.ndata = ivector(1,wpx.nsamples);
  wpx.ndatar = ivector(1,wpx.nsamples);
  wpx.ndatab = ivector(1,wpx.nsamples);
  wpx.ngal = dvector(1,wpx.nsamples);

  
  set_lensing_mass_bins();

  // loop over BLUE and RED
  for(ii=1;ii<=2;++ii)
    {
      if(ii==1)sprintf(color,"blue");
      if(ii==2)sprintf(color,"red");

      // read in the wp data
      for(i=istart;i<=wpx.nsamples;++i)
	{
	  // Tinker's laptop
	  sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z%d_%s_sm%d.fits.txt",
		  wpx.iz,color,i-1);
	  if(JPL_FLAG)
	    sprintf(fname,"./DATA/z%d_%s_sm%d.fits.txt",wpx.iz,color,i-1);

	  fp = openfile(fname);
	  n = filesize(fp);
	  wpx.ndata[i] = n;
	  if(ii==1)
	    {
	      wpx.rdata[i] = dvector(1,nlim);
	      wpx.xdata[i] = dvector(1,nlim);
	      wpx.edata[i] = dvector(1,nlim);
	      wpx.covar[i] = dmatrix(1,nlim,1,nlim);
	      wpx.model[i] = dvector(1,nlim);
	    }

	  for(j=1;j<=n;++j)
	    {
	      fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		     &x1,&x1,&x1,&x1,&wpx.rdata[i][j],&x1,&wpx.xdata[i][j],&x1,&wpx.edata[i][j]);
	      fgets(aa,1000,fp);
	    }
	  fprintf(stderr,"Read %d lines from [%s]\n",n,fname);
	  fclose(fp);

	  /* let's take these data and put them in comoving h-inverse Mpc
	   * r_ph = r_co*a = r_co/(1+z)
	   */
	  for(j=1;j<=n;++j)
	    {
	      wpx.rdata[i][j] = wpx.rdata[i][j]/1000.0*HUBBLE*(1+REDSHIFT);
	      //printf("%d %d %e %e %e\n",i,j,wpx.rdata[i][j],wpx.xdata[i][j],wpx.edata[i][j]);
	    }

	  if(LENSING_COVAR)
	    {
	      // Now read in the covariance matrix
	      // Tinker's laptop
	      sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/Covar_%s/z%d_%s_sm%d.covar",
		      color,wpx.iz,color,i-1);
	      if(JPL_FLAG)
		sprintf(fname,"./DATA/z%d_%s_sm%d.covar",wpx.iz,color,i-1);
	      
	      fp = openfile(fname);
	      if(filesize(fp)!=n*n)
		{
		  fprintf(stderr,"FILESIZE mismatch: cij and data not match in SMF: %d %d\n",n*n,filesize(fp));
		  exit(0);
		}
	      covar = dmatrix(1,nlim,1,nlim);
	      for(j=1;j<=n;++j)
		for(k=1;k<=n;++k)
		  fscanf(fp,"%d %d %lf",&i1,&i1,&covar[j][k]);
	      fclose(fp);	
	      fprintf(stderr,"Read [%d] lines from [%s]\n",n*n,fname);
	      
	      //add the statistical uncertainties to the diagonals
	      for(j=1;j<=n;++j)
		covar[j][j] += (wpx.edata[i][j]*wpx.edata[i][j]);


	      printf("INVERTING COVARIANCE MATRIX\n");
	      tmp=dmatrix(1,n,1,1);
	      tmp2=dmatrix(1,n,1,n);
	      for(j=1;j<=n;++j)
		for(k=1;k<=n;++k)
		  tmp2[j][k]=covar[j][k];
	      gaussj(tmp2,n,tmp,1);
	      for(j=1;j<=n;++j)
		for(k=1;k<=n;++k)
		  wpx.covar[i][j][k]=tmp2[j][k];
	      free_dmatrix(tmp,1,n,1,1);
	      free_dmatrix(tmp2,1,n,1,n);
	      free_dmatrix(covar,1,n,1,n);
	    }
	  // put this in the right place
	  if(ii==1)
	    {
	      wpx.ndatab[i] = n;
	      wpx.rdatab[i] = dvector(1,nlim);
	      wpx.xdatab[i] = dvector(1,nlim);
	      wpx.edatab[i] = dvector(1,nlim);
	      wpx.covarb[i] = dmatrix(1,nlim,1,nlim);
	      wpx.two_halo_blue[i] = dvector(1,nlim);

	      for(j=1;j<=n;++j)
		{
		  wpx.rdatab[i][j] = wpx.rdata[i][j];
		  wpx.xdatab[i][j] = wpx.xdata[i][j];
		  wpx.edatab[i][j] = wpx.edata[i][j];
		  for(k=1;k<=n;++k)
		    wpx.covarb[i][j][k] = wpx.covar[i][j][k];
		}
	    }
	  if(ii==2)
	    {
	      wpx.ndatar[i] = n;
	      wpx.rdatar[i] = dvector(1,nlim);
	      wpx.xdatar[i] = dvector(1,nlim);
	      wpx.edatar[i] = dvector(1,nlim);
	      wpx.covarr[i] = dmatrix(1,nlim,1,nlim);
	      wpx.two_halo_red[i] = dvector(1,nlim);

	      for(j=1;j<=n;++j)
		{
		  wpx.rdatar[i][j] = wpx.rdata[i][j];
		  wpx.xdatar[i][j] = wpx.xdata[i][j];
		  wpx.edatar[i][j] = wpx.edata[i][j];
		  for(k=1;k<=n;++k)
		    wpx.covarr[i][j][k] = wpx.covar[i][j][k];
		}
	    }
	}
    }
}

/****
 */
void set_lensing_mass_bins()
{
  int i;
  if(wpx.iz==1)
    {
      // bins for the BLUE SAMPLE
      wpx.mstar_wplo[1][1] = 11.12;
      wpx.mstar_wplo[1][2] = 10.89;
      wpx.mstar_wplo[1][3] = 10.64;
      wpx.mstar_wplo[1][4] = 10.3;
      wpx.mstar_wplo[1][5] = 9.82;
      wpx.mstar_wplo[1][6] = 9.2;
      wpx.mstar_wplo[1][7] = 8.7;
      
      wpx.mstar_wphi[1][1] = 12.0;
      wpx.mstar_wphi[1][2] = 11.12;
      wpx.mstar_wphi[1][3] = 10.89;
      wpx.mstar_wphi[1][4] = 10.64;
      wpx.mstar_wphi[1][5] = 10.3;
      wpx.mstar_wphi[1][6] = 9.82;
      wpx.mstar_wphi[1][7] = 9.2;
      
      // bins for the RED SAMPLE (same as blue sample)
      wpx.mstar_wplo[2][1] = 11.12;
      wpx.mstar_wplo[2][2] = 10.89;
      wpx.mstar_wplo[2][3] = 10.64;
      wpx.mstar_wplo[2][4] = 10.3;
      wpx.mstar_wplo[2][5] = 9.82;
      wpx.mstar_wplo[2][6] = 9.2;
      wpx.mstar_wplo[2][7] = 8.7;
      
      wpx.mstar_wphi[2][1] = 12.0;
      wpx.mstar_wphi[2][2] = 11.12;
      wpx.mstar_wphi[2][3] = 10.89;
      wpx.mstar_wphi[2][4] = 10.64;
      wpx.mstar_wphi[2][5] = 10.3;
      wpx.mstar_wphi[2][6] = 9.82;
      wpx.mstar_wphi[2][7] = 9.2;
    }

  if(wpx.iz==2)
    {
      // bins for the BLUE SAMPLE
      wpx.mstar_wplo[1][1] = 11.29;
      wpx.mstar_wplo[1][2] = 11.05;
      wpx.mstar_wplo[1][3] = 10.88;
      wpx.mstar_wplo[1][4] = 10.65;
      wpx.mstar_wplo[1][5] = 10.3;
      wpx.mstar_wplo[1][6] = 9.8;
      wpx.mstar_wplo[1][7] = 9.3;
      
      wpx.mstar_wphi[1][1] = 12.0;
      wpx.mstar_wphi[1][2] = 11.29;
      wpx.mstar_wphi[1][3] = 11.05;
      wpx.mstar_wphi[1][4] = 10.88;
      wpx.mstar_wphi[1][5] = 10.65;
      wpx.mstar_wphi[1][6] = 10.3;
      wpx.mstar_wphi[1][7] = 9.8;

      for(i=1;i<=7;++i)
	wpx.mstar_wplo[2][i] = wpx.mstar_wplo[1][i];
      for(i=1;i<=7;++i)
	wpx.mstar_wphi[2][i] = wpx.mstar_wphi[1][i];
    }

  if(wpx.iz==3)
    {
      // bins for the BLUE SAMPLE
      wpx.mstar_wplo[1][1] = 11.35;
      wpx.mstar_wplo[1][2] = 11.16;
      wpx.mstar_wplo[1][3] = 10.97;
      wpx.mstar_wplo[1][4] = 10.74;
      wpx.mstar_wplo[1][5] = 10.39;
      wpx.mstar_wplo[1][6] = 9.8;
      
      wpx.mstar_wphi[1][1] = 12.0;
      wpx.mstar_wphi[1][2] = 11.35;
      wpx.mstar_wphi[1][3] = 11.16;
      wpx.mstar_wphi[1][4] = 10.97;
      wpx.mstar_wphi[1][5] = 10.74;
      wpx.mstar_wphi[1][6] = 10.39;

      for(i=1;i<=6;++i)
	wpx.mstar_wplo[2][i] = wpx.mstar_wplo[1][i];
      for(i=1;i<=6;++i)
	wpx.mstar_wphi[2][i] = wpx.mstar_wphi[1][i];
    }


}
