#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* Calculate and output the galaxy-galaxy clustering signal 
 * in the user-defined stellar mass bins.
 */
void shmr_clustering()
{
  int i,j;
  double mlo, mhi, x1, x2, x3, x4, r, dr;
  FILE *fp;
  char fname[1000];

  if(COLOR)return shmr_color_clustering();

  wpl.reset_inversion = 1;
  for(i=1;i<=wpl.nsamples;++i)
    {
      sprintf(fname,"%s.clustering_bin%d",Task.root_filename,i);
      fp = fopen(fname,"w");
      // make a header:
      fprintf(fp,"# %e %e\n",wpl.mstar_wplo[0][i],wpl.mstar_wphi[0][i]);
      fflush(fp);
      mlo = pow(10.0,wpl.mstar_wplo[0][i]);
      mhi = pow(10.0,wpl.mstar_wphi[0][i]);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);

      GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      fprintf(stdout,"%e %e\n",GALAXY_DENSITY,GALAXY_BIAS);
      RESET_FLAG_1H = RESET_FLAG_2H = 1;
      RESET_KAISER++;

      dr=(log(70.0)-log(0.01))/49.0;
      for(j=0;j<50;++j)
	{
	  r = exp(j*dr+log(0.01));
	  x1 = one_halo_real_space(r);
	  x2 = two_halo_real_space(r);
	  x3 = projected_xi(r);
	  x4 = projected_xi_1halo(r);
	  fprintf(fp,"%f %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4);
	  fflush(fp);
	}
      fclose(fp);      
    }
}

void shmr_color_clustering(void)
{  
  int i,j;
  double mlo, mhi, x1, x2, x3, x4, r, dr;
  FILE *fp;
  char fname[1000];

  BLUE_FLAG = 1;
  for(i=1;i<=wpl.nsamples;++i)
    {
      wpl.reset_inversion = 1;
      RESET_FLAG_1H = 1;
      RESET_FLAG_2H = 1;
      sprintf(fname,"%s.blue_clustering_bin%d",Task.root_filename,i);
      fp = fopen(fname,"w");
      // make a header:
      fprintf(fp,"# %e %e\n",wpl.mstar_wplo[0][i],wpl.mstar_wphi[0][i]);
      fflush(fp);
      mlo = pow(10.0,wpl.mstar_wplo[0][i]);
      mhi = pow(10.0,wpl.mstar_wphi[0][i]);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);
      dr=(log(70.0)-log(0.01))/49.0;
      for(j=0;j<50;++j)
	{
	  r = exp(j*dr+log(0.01));
	  x1 = one_halo_real_space(r);
	  x2 = two_halo_real_space(r);
	  x3 = projected_xi(r);
	  x4 = 0;//projected_xi_1halo(r);
	  fprintf(fp,"%f %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4);
	  fflush(fp);
	}
      fclose(fp);      
    }

  BLUE_FLAG = 0;
  for(i=1;i<=wpl.nsamples;++i)
    {
      wpl.reset_inversion = 1;
      RESET_FLAG_1H = 1;
      RESET_FLAG_2H = 1;
      sprintf(fname,"%s.red_clustering_bin%d",Task.root_filename,i);
      fp = fopen(fname,"w");
      // make a header:
      fprintf(fp,"# %e %e\n",wpl.mstar_wplo[0][i],wpl.mstar_wphi[0][i]);
      fflush(fp);
      mlo = pow(10.0,wpl.mstar_wplo[0][i]);
      mhi = pow(10.0,wpl.mstar_wphi[0][i]);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);
      dr=(log(70.0)-log(0.01))/49.0;
      for(j=0;j<50;++j)
	{
	  r = exp(j*dr+log(0.01));
	  x1 = one_halo_real_space(r);
	  x2 = two_halo_real_space(r);
	  x3 = projected_xi(r);
	  x4 = projected_xi_1halo(r);
	  fprintf(fp,"%f %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4);
	  fflush(fp);
	}
      fclose(fp);      
    }

  return;
}

/* Read in the setllar mass bins from the
 * user-supplied file listed in the hod.bat file.
 */
void input_stellar_mass_bins()
{
  FILE *fp;
  int i,j,n;

  if(!Task.clustering && !Task.populate_sim)goto INPUT_LENSING_FILES;
  if(!Files.UseStellarMassBinsClustering)
    {
      fprintf(stdout,"No StellarMassBins file specified for clustering, using defaults:\nlog(M0) = 9, with 0.5 dex increments up to 12.0\n");
      wpl.nsamples = (12-9)/0.5;
      for(i=1;i<=wpl.nsamples;++i)
	{
	  wpl.mstar_wplo[0][i] = 9.0 + 0.5*i;
	  wpl.mstar_wphi[0][i] = 9.5 + 0.5*i;
	}
    }
  else
    {
      fprintf(stdout,"Opening file: [%s]\n",Files.StellarMassBinsClustering);
      fp = openfile(Files.StellarMassBinsClustering);
      wpl.nsamples = filesize(fp);
      
      printf("Reading %d stellar mass bins.\n",wpl.nsamples);
      for(i=1;i<=wpl.nsamples;++i)
	fscanf(fp,"%lf %lf",&wpl.mstar_wplo[0][i], &wpl.mstar_wphi[0][i]);
      fclose(fp);
    }

 INPUT_LENSING_FILES:
  if(!Task.lensing)return;
  if(!Files.UseStellarMassBinsLensing)
    {
      fprintf(stdout,"No StellarMassBins file specified for lensing, using defaults:\nlog(M0) = 9, with 0.5 dex increments up to 12.0\n");
      wpx.nsamples = (12-9)/0.5;
      for(i=1;i<=wpx.nsamples;++i)
	{
	  wpx.mstar_wplo[0][i] = 9.0 + 0.5*i;
	  wpx.mstar_wphi[0][i] = 9.5 + 0.5*i;
	}
    }
  else
    {
      fprintf(stdout,"Opening file: [%s]\n",Files.StellarMassBinsLensing);
      fp = openfile(Files.StellarMassBinsLensing);
      wpx.nsamples = filesize(fp);
      
      printf("Reading %d stellar mass bins.\n",wpx.nsamples);
      for(i=1;i<=wpx.nsamples;++i)
	fscanf(fp,"%lf %lf",&wpx.mstar_wplo[0][i], &wpx.mstar_wphi[0][i]);
      fclose(fp);
    }

}

void shmr_lensing()
{
  int i,j;
  double mlo, mhi, x1, x2, x3, x4, r, dr;
  FILE *fp;
  char fname[1000];

  if(COLOR)return shmr_color_lensing();

  for(i=1;i<=wpx.nsamples;++i)
    {
      wpx.reset_inversion = 1;
      RESET_FLAG_XGM1 = 1;
      RESET_FLAG_2H = 1;
      RESET_FLAG_DS = 1;
      sprintf(fname,"%s.lensing_bin%d",Task.root_filename,i);
      fp = fopen(fname,"w");
      // make a header:
      fprintf(fp,"# %e %e\n",wpx.mstar_wplo[0][i],wpx.mstar_wphi[0][i]);
      fflush(fp);
      mlo = pow(10.0,wpx.mstar_wplo[0][i]);
      mhi = pow(10.0,wpx.mstar_wphi[0][i]);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);
      dr=(log(50.0)-log(0.01))/49.0;
      for(j=0;j<50;++j)
	{	  
	  r = exp(j*dr+log(0.01));
	  x3 = delta_sigma(r);
	  x1 = one_halo_galaxy_matter(r);
	  x2 = two_halo_real_space(r);
	  fprintf(fp,"%f %e %e %e %e\n",r,x1,x2,x1+x2,x3);
	  fflush(fp);
	}
      fclose(fp);      
    }
}

void shmr_color_lensing()
{
}

