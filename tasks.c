#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

/***********************************************************************
 * This is the main function for controlling the activities of the code.
 * If you ask for any sort of minimization or MCMC analysis, the code is
 * redirected to another routine, but once finished will complete the rest
 * of the tasks checked off (so if you want to minimize and then print out
 * related statistics, the code will do that automatically if you ask it).
 *
 *
 * Each task generates an output file using the "root_filename" in the hod.bat
 * file and tacking on an extension at the end:
 *
 *  [filename].r_space  --> real space correlation function
 *  [filename].r_half   --> r_xi/2 statistic (z-space)
 *  [filename].qp_ratio --> multipole diagnostics of z-space xi
 *  [filename].linlin   --> 2-d z-space xi with linear bins
 *  [filename].loglog   --> 2-d z-space xi with logarithmic bins
 *  [filename].HOD      --> prints out the mean N(M)
 */


void tasks(int argc, char **argv)
{
  int i,j,nsize,ix,iy,nr,i1,i2,n,ibuf;
  float x1,x2,x3,x4,x5,x6,err,npairs;
  double r,rs,rp,dx1,dx2,dx3,dx4,dx5,**xi2d,*xx2d,*yy2d,**xi2d_data,**xi2d_kaiser,
    xi0_m,xi2_m,xi0_k,xi2_k,xi0_d,xi2_d,xi_r,rlo,rhi,rr[50],rhalf[50],dr,
    rphi,rshi,rslo,rplo,dlogm,m,sig,mhi,mlo,ms,fsat,fsat_prev;
  FILE *fp,*fp2,*fp1;
  char fname[100];

  
  /* If the bias if to be fit for from the clustering_filenames.dat
   * files...
   */
  if(Task.fit_for_bias)
    fit_for_bias();

  /* This is for chi^2 minimization of data for the projected correlation
   * function.
   */
  if(Task.wp_minimize)
    wp_minimization(argv[1]);

  /* This is for chi^2 minimization of the combined data for stellar-to-halo mass relation
   */
  if(Task.shmr_minimize)
    shmr_minimization();

  /* This is for Monte-Carlo Markov Chain exploration of the posterior
   * distribution of the parameters, either real-space or redshift-space,
   * depending on what MCMC is set to.
   */
  if(Task.MCMC==1 && SHMR_FLAG==0) 
    mcmc_minimization();

  if(Task.MCMC==1 && SHMR_FLAG==1) 
    mcmc_shmr();

  if(Task.MCMC==3 && SHMR_FLAG==1) 
    chain_postprocess();

  /* This is to output the shape of the mean occupation function and the 
   * scatter about the mean. 
   * File columns for [root].HOD are:
   *  1 - halo mass (M_sol/h)
   *  2 - N_cen(M)
   *  3 - N_sat(M)
   *  4 - N_tot(M)
   *  5 - <N(N-1)> 
   *  6 - ratio of scatter to Poisson scatter (often called "alpha")
   * NB! if SHMR_FLAG set, each stellar mass bin is outputted in its own
   * file, [root].HOD_bin[i]
   */
  if(Task.HOD)
    {
      if(!SHMR_FLAG) { 
	sprintf(fname,"%s.HOD",Task.root_filename);
	fp=fopen(fname,"w");
	dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
	for(i=1;i<=100;++i)
	  {
	    m=exp((i-1)*dlogm)*HOD.M_low;
	    sig = N_sat(m)*N_sat(m) + N_cen(m)*(1-N_cen(m));
	    fprintf(fp,"%e %e %e %e %e %e\n",
		    m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)));
	  }
	fclose(fp);
      }  else {
	for(i=1;i<=wpl.nsamples;++i)
	  {
	    mlo = pow(10.0,wpl.mstar_wplo[0][i]);
	    mhi = pow(10.0,wpl.mstar_wphi[0][i]);
	    set_up_hod_for_shmr(mlo,mhi,wpl.a);
	    sprintf(fname,"%s.HOD_bin%d",Task.root_filename,i);
	    fp = fopen(fname,"w");
	    // make a header:
	    fprintf(fp,"# %e %e\n",wpl.mstar_wplo[0][i],wpl.mstar_wphi[0][i]);
	    fprintf(fp,"# %e %e %e %e %e\n",HOD.M_min, HOD.M1, HOD.M_cut, HOD.alpha, wpl.a[6]);
	    dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
	    for(j=1;j<=100;++j)
	      {
		m=exp((j-1)*dlogm)*HOD.M_low;
		sig = N_sat(m)*N_sat(m) + N_cen(m)*(1-N_cen(m));
		fprintf(fp,"%e %e %e %e %e %e\n",
			m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)));
	      }
	    fclose(fp);
	  }
      }
      
    }

  if(Task.lensing)
    {
      fprintf(stderr,"\n\nCALCULATING GALAXY-MATTER CROSS CORRELATION FUNCTION.\n");
      fprintf(stderr,    "----------------------------------------------------\n\n");
      shmr_lensing();
    }

  /* This sets in motion the calculation and tabulation of the real-space
   * correlation function and outputs it to a file.
   * File columns for [root].clustering are:
   *  1 - radius (Mpc/h)
   *  2 - one-halo term (real-space)
   *  3 - two-halo term (real-space)
   *  4 - full xi(r)
   *  5 - projected correlation function (1st column is now r_p)
   *  6 - w_p of the one-halo term only
   * NB! - if SHMR_FLAG set, output comes from shmr_clustering.c
   */
  if(Task.clustering)
    {
      fprintf(stderr,"\n\nCALCULATING REAL-SPACE CORRELATION FUNCTION.\n");
      fprintf(stderr,    "--------------------------------------------\n\n");
      if(SHMR_FLAG)shmr_clustering();
      else {
	sprintf(fname,"%s.clustering",Task.root_filename);
	fp=fopen(fname,"w");
	dr=(log(70.0)-log(0.01))/49.0;
	for(i=0;i<50;++i)
	  {
	    r=exp(i*dr+log(0.01));
	    x1=one_halo_real_space(r);
	    x2=two_halo_real_space(r);
	    x3=projected_xi(r);
	    x4 = projected_xi_1halo(r);
	    x5 = projected_xi_rspace(r);
	    fprintf(fp,"%f %e %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4,x5);
	    fflush(fp);
	  }
	fclose(fp);
      }
    }

  /* This outputs the SHMR that is specified in the hod.bat file
   * Filename: [root].SHMR
   * File columns are:
   *  1 - Stellar Mass [Msol]
   *  2 - Halo Mass [Msol/h]
   *  3 - satellite fraction
   *  4 - Stellar mass function [(Mpc/h)^-3 dex^-1]
   *  5 - Galaxy bias
   * Apologies for the weird units.
   */
  if(Task.SHMR && SHMR_FLAG && !COLOR)
    {
      fprintf(stderr,"\n\nCALCULATING STELLAR-TO-HALO MASS RELATION.\n");
      fprintf(stderr,    "------------------------------------------\n\n");
      sprintf(fname,"%s.SHMR",Task.root_filename);
      fp = fopen(fname,"w");
      mlo = 1.0E9;
      mhi = 5.0e12;
      dlogm = log(mhi/mlo)/49;      
      for(i=0;i<50;++i)
	{
	  ms = exp(i*dlogm)*mlo;
	  set_up_hod_for_shmr(mlo*exp((i-0.5)*dlogm),mlo*exp((i+0.5)*dlogm),wpl.a);
	  fsat = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	  fprintf(fp,"%e %e %e %e %e %e %e\n",ms,ms_to_mhalo(ms,wpl.a),fsat,GALAXY_DENSITY/dlogm*log(10),GALAXY_BIAS,HOD.M1/HOD.M_min,HOD.M_min);
	}
      fclose(fp);
    }

  if(Task.SHMR && SHMR_FLAG && COLOR)
    {
      fsat_prev = 1;
      fprintf(stderr,"\n\nCALCULATING STELLAR-TO-HALO MASS RELATION.\n");
      fprintf(stderr,    "------------------------------------------\n\n");
      BLUE_FLAG = 1;
      sprintf(fname,"%s.SHMR_blue",Task.root_filename);
      fp = fopen(fname,"w");
      mlo = 1.0E9;
      mhi = 1.0e12;
      dlogm = log(mhi/mlo)/49;      
      for(i=0;i<50;++i)
	{
	  ms = exp(i*dlogm)*mlo;
	  set_up_hod_for_shmr(mlo*exp((i-0.5)*dlogm),mlo*exp((i+0.5)*dlogm),wpl.a);
	  fsat = 0;
	  if(fsat_prev>1.0E-6) {
	    fsat = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	    fsat_prev = fsat;
	  }
	  fprintf(fp,"%e %e %e %e %e\n",ms,ms_to_mhalo(ms,wpl.a),fsat,GALAXY_DENSITY/dlogm*log(10),GALAXY_BIAS);
	  fflush(fp);
	}
      fclose(fp);

      BLUE_FLAG = 0;
      fsat_prev = 1;
      ibuf = SHMR_PARAMS;
      wpl.reset_inversion = 1;
      sprintf(fname,"%s.SHMR_red",Task.root_filename);
      fp = fopen(fname,"w");
      mlo = 1.0E9;
      mhi = 1.0e12;
      dlogm = log(mhi/mlo)/49;      
      for(i=0;i<50;++i)
	{
	  ms = exp(i*dlogm)*mlo;
	  set_up_hod_for_shmr(mlo*exp((i-0.5)*dlogm),mlo*exp((i+0.5)*dlogm),wpl.a);
	  fsat = 0;
	  if(fsat_prev>1.0E-6) {
	    fsat = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	    fsat_prev = fsat;
	  }
	  fprintf(fp,"%e %e %e %e %e\n",ms,ms_to_mhalo(ms,&wpl.a[ibuf]),fsat,GALAXY_DENSITY/dlogm*log(10),GALAXY_BIAS);
	}
      fclose(fp);
    }

  /* This takes a halofile from a simulation and populates the halos
   * with galaxies according the HOD specified in the batch file.
   */ 
  if(Task.populate_sim==1 && SHMR_FLAG)
    populate_simulation();
  if(Task.populate_sim==1 && !SHMR_FLAG)
    populate_simulation_hod();

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear P(k) calculated using Smith et al.
   *
   * Format of file [root].matter_pk
   *   1 - k [h/Mpc]
   *   2 - linear P(k) [Mpc/h]^3
   *   3 - non-linear P(k) [Mpc/h]^3
   */
  if(Task.matter_pk)
    output_matter_power_spectrum();

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear xi(r) is Fourier transform of Smith et al (above)
   *
   * Format of file [root].matter_pk
   *   1 - r [Mpc/h]
   *   2 - linear xi(r)
   *   3 - non-linear xi(r)
   */
  if(Task.matter_xi)
    output_matter_correlation_function();

  /* Output the matter variance as a function of scale.
   *
   * Format of file [root].sigma_r 
   *  1 - r [Mpc/h]
   *  2 - sigma(r) [linear]
   *  3 - sigma(r) [non-linear, using Smith]
   *  4 - mass [M_sol/h] mass = (4/3)*PI*r^3*rho_bar
   */
  if(Task.sigma_r)
    output_matter_variance();

  /* Output the halo concentrations
   * 
   * Format of file [root].civr
   *  1 - mass [Msol/h] --> mass specified by DELTA_HALO (input file).
   *  2 - halo concentration.
   */
  if(Task.cvir)
    output_halo_concentrations();

  /* Output halo abundance and bias
   * 
   * Format of file [root].halostats
   *  1 - mass [Msol/h] --> mass specified by DELTA_HALO (input file).
   *  2 - halo abundance, dn/dM, [1/(Msol/h)/(Mpc/h)^3].
   *  3 - halo bias (no units)
   */
  if(Task.halostats)
    output_halo_mass_function();

}
