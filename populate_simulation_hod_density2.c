#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

#define NBINLOOKUP 10000

/* This function reads in a halo file and creates a mock galaxy distribution 
 * by populating the dark matter halos with galaxies with the currently specified
 * HOD functions.
 *
 * The format of the halo file needs to be in ASCII:
 * col 1 - halo ID (not used)
 * col 2 - number of particles in halo
 * cols 3-5 - x,y,z position of halo, in Mpc/h
 * col 6 - space velocity of halo (not used)
 * cols 7-9 - vx,vy,vz velocity of halo, in km/s
 *
 * The values of RESOLUTION, BOX_SIZE, OMEGA_M from the batch file are used
 * to convert particle number into mass, and to box-wrap galaxies.
 * If you want another file format, by all means edit.
 *
 * Output: mock galaxies are printed to [root].mock
 * in format: x,y,z [Mpc/h] vz,vy,vz [km/s]
 *
 * Output: HOD calculated from populated halos (for cross-checking with true HOD)
 * [root].binned_HOD in format: bin id, log10(mass), <N_gal>, N_halo
 * 
 * Satellite positions are chosen from a random sampling of the NFW profile
 * for the mass of the halo being populated. If CVIR_FAC is specified, then that value
 * will be used to adjust the NFW profile. Satellite velocities are chosen
 * from a Gaussin distribution with width = to virial dispersion of halo (plus
 * halo center-of-mass).
 *
 * NB-- If you are using the code to calculate M_min from a desired space density,
 * be sure that the linear matter power spectrum is the same as that used in the
 * simulation, or else the space densities won't match. [Mass function is used for this
 * purpose].
 *
 */
double *g21_rad, *g21_xi;
int g21_nrad;
extern FILE *FP1;

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);
void mmin_from_file(void);


/* Internal function
 */
void calc_nbody_two_halo(float **gal, int *id, int ngal);
void M_min_delta(double m, double delta, double norm);
double find_norm(double norm);

double box_n_glob_pop;
double slope_glob_pop=-0.10, delta_crit_glob_pop=1.6;
float *mass_glob_pop, *delta_glob_pop;
int nhalos_glob_pop;

void M_min_delta_challenge(double m, double delta, double norm)
{
  double factor;

  factor = norm + slope_glob_pop * (delta-1.6);
  if (factor < 0.01)
    HOD.M_min = HOD.M_min_0 * 0.01;
  else
    HOD.M_min = HOD.M_min_0 * factor;
}

/* changing this function for the Fisher emulator calculations.
 */
void M_min_delta(double m, double delta, double norm)
{
  double factor;
  
  factor = norm + (1 + HOD.f_env*erf((delta-HOD.delta_env)/HOD.sigma_env));
  if (factor < 0.01) {
    HOD.M_low = HOD.M_low_0 * 0.01;
    HOD.M_cut = HOD.M_cut_0 * 0.01;
    HOD.M_min = HOD.M_min_0 * 0.01; }
  else {
    HOD.M_low = HOD.M_low_0 * factor;
    HOD.M_cut = HOD.M_cut_0 * factor;
    HOD.M_min = HOD.M_min_0 * factor; }
}

double find_norm(double norm)
{
  double mass, ncen, ngals=0, nsat;
  int i;

  for (i = 0; i < nhalos_glob_pop; ++i)
    {
      mass = mass_glob_pop[i];
      M_min_delta(mass, delta_glob_pop[i], norm);
      if(mass > HOD.M_max) continue;
      if(mass < HOD.M_low) continue;

      //printf("%e %e %e %e %e\n",HOD.M_min_0, HOD.M_min,delta_glob_pop[i],mass,N_cen(mass));
      ncen=N_cen(mass);
      ngals += ncen;
      nsat = N_sat(mass);
      ngals += nsat;
      if(nsat<0) { fmuh(mass); fmuh(ngals); fmuh(nsat); fmuh(ncen); }
      if(ncen<0) { fmuh(mass); fmuh(ngals); fmuh(nsat); fmuh(ncen); }
      if(isnan(ngals)){ fmuh(mass); exit(0); }
    }
  ngals /= (BOX_SIZE*BOX_SIZE*BOX_SIZE);

  fprintf(stdout,"\nNORMALISATION FACTOR:\t%f\nNUMBER DENSITY OF GALAXIES:\t%e\ndifference:\t%e\n", norm, ngals, 
	  (ngals - box_n_glob_pop)/box_n_glob_pop);
  fflush(stdout);
  //exit(0);
  return ((ngals - box_n_glob_pop)/box_n_glob_pop);
  
}

void populate_simulation_hod_density()
{
  FILE *fp, *fdelta, *fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt, *fgal, *fp1;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  int nhalos, haloid;
  float *xhs, *yhs, *zhs, *vxhs, *vyhs, *vzhs, rootnorm;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv, ngals;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax, omega_z, vfac;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0;
  int HOD_ENV=1, ii;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  int SO_FILE = 1,
    JEANS_DISPERSION = 0;

  //  if(RESOLUTION > 1.6)
  // WARREN_MASS_CORRECTION = 1;

  TRACK_GALAXIES=0;
  ngals=0.0;

  galarr = matrix(1,10000000,0,5);
  haloarr = matrix(1,10000000,0,6);
  galid = ivector(1,10000000);

  fp=openfile(Files.HaloFile);
  //fdelta = openfile("halo_density.0029.dat");
  //fdelta = openfile("/data1/tinker/HALO_CATALOGS/BigMULTIDARK_Planck/halo_density_2e12.dat");
  sprintf(aa,"%s.delta",Files.HaloFile);
  fdelta = openfile(aa);

  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");
  //sprintf(aa,"%s.mock_halo",Task.root_filename);      
  //fpt = fopen(aa,"w");

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(i=0;i<1000;++i)
    galcnt[i]=0;

  // NOT SO FAST! let's get the hod params from the box itself
  //set_HOD_params();
  fp1 = openfile("omega_value.dat");
  fscanf(fp1,"%lf",&OMEGA_M);
  fclose(fp1);
  fprintf(stdout,"NEW OMEGA_M: %f\n",OMEGA_M);

  // read in the whole HOD from a file
  fp1 = openfile(ARGV[3]);
  fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	 &GAMMA, &HOD.M1, &HOD.alpha, &HOD.M_cut, &HOD.sigma_logM, 
	 &CVIR_FAC, &VBIAS, &VBIAS_C, &HOD.f_env, &HOD.delta_env, &HOD.sigma_env);
  fprintf(stdout,"NEW HOD: %e %e %e %e %e %e %e %e\n",GAMMA, HOD.M1,HOD.alpha, HOD.M_cut, HOD.sigma_logM, 
	 HOD.f_env, HOD.delta_env, HOD.sigma_env);
  fprintf(stdout,"NEW HOD: %e %e %e\n",CVIR_FAC,VBIAS,VBIAS_C);
  OMEGA_Z = OMEGA_M*1.55*1.55*1.55/(OMEGA_M*1.55*1.55*1.55 + (1-OMEGA_M));
  vfac = pow(OMEGA_Z,GAMMA)/pow(OMEGA_Z,0.55);      
  fprintf(stdout,"VFAC: %e\n", vfac);

  /*
  // read in the density params
  HOD.f_env = atof(ARGV[3]);
  HOD.sigma_env = atof(ARGV[4]);
  HOD.delta_env = atof(ARGV[5]);
  fprintf(stdout,"HOD_env: %f %f %f\n",HOD.f_env, HOD.sigma_env, HOD.delta_env);
  // if more params, read them in
  if(ARGC>6)
    {
      HOD.sigma_logM = atof(ARGV[6]);
      fprintf(stdout,"HOD.slogm = %f\n",HOD.sigma_logM);
    }
  */

  FP1 = fp;
  //ZHONGXU!!
  mmin_from_file(); // only use this if you are GALAXY_DENSITY SET
  // if you want to set M_min and determine GALAXY_DENSITY, you need to write new code
  // that takes the input HOD parameters (the global ones)
  // and finds GALAXY_DENSITY from that. should be easier than finding M_min itself.
  // GALAXY_DENSITY = function_call_find_density();
  /*
  Ntot = 0;
  for(i=0;i<Nhalos;++i)
    { 
      fscanf(fphalofile, "%lf", &halomass);
      Ntot += N_cen(halomass) + N_sat(halomass);
    }
  GALAXY_DENSITY = Ntot/pow(BOX_SIZE,3.0);
  */   

  HOD.M_low_0 = HOD.M_low;
  HOD.M_min_0 = HOD.M_min;
  HOD.M_cut_0 = HOD.M_cut;
  mlo = HOD.M_low;
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  printf("BOX_SIZE %f\n",BOX_SIZE);
  fflush(stdout);

  haloid = 0;

  
  /* DENSITY CALC!!
   * calculate the number density now by looping over the halos and 
   * imparting some density dependence on the hod
   */
  // step 1: read the halo file into memory (completely)
  // step 2: read the density file into memory
  // loop over halo mass.

  // HOD functions
  // N_cen(mass)
  // N_sat(mass)

  // for i in all halos
  // ntot += N_cen(masss[i]) + N_sat(mas[i])

  nhalos_glob_pop = filesize(fp);
  nhalos = nhalos_glob_pop;
  xhs = vector(0, nhalos - 1);
  yhs = vector(0, nhalos - 1);
  zhs = vector(0, nhalos - 1);
  vxhs = vector(0, nhalos - 1);
  vyhs = vector(0, nhalos - 1);
  vzhs = vector(0, nhalos - 1);

  mass_glob_pop = vector(0, nhalos - 1);
  delta_glob_pop = vector(0, nhalos - 1);
  
  if(filesize(fdelta)!= nhalos) 
    {
      fprintf(stderr,"ERROR: filesize mismatch of halos with density!\n");
      exit(0);
    }
  
  if(MCMC==3)
    {
      omega_z = 0.30711*pow(1+REDSHIFT,3.0)/(0.30711*pow(1+REDSHIFT,3.0) + (1-0.30711));
      printf("populate_simulation_hod_density> mfac=%e vfac=%e\n",(OMEGA_TEMP/0.30711),vfac);
    }
  for(i = 0; i < nhalos; ++i)
    {
      // TEMP! this is for the PM emulator with FOF halos
      // ZHONGXU-- make sure read statement is compatible with the halo catalogs you are using!!!
      fscanf(fp,"%d %d %f %f %f %f %f %f %f",
	     &j, &imass,&xhs[i], &yhs[i], &zhs[i], &x1, &vxhs[i], &vyhs[i], &vzhs[i]);
      mass_glob_pop[i] = imass*pow(RESOLUTION,3)*OMEGA_M*RHO_CRIT;

      
      vxhs[i] *= vfac;
      vyhs[i] *= vfac;
      vzhs[i] *= vfac;

      //fscanf(fp,"%f %f %f %f %f %f %f %f %f %f",
      //     &x1, &mass_glob_pop[i],&x1,&x1,&xhs[i], &yhs[i], &zhs[i], &vxhs[i], &vyhs[i], &vzhs[i]);
      //fscanf(fp,"%f %f %f %f %f %f %f %f",
      //       &xhs[i], &yhs[i], &zhs[i], &vxhs[i], &vyhs[i], &vzhs[i], &x1, &mass_glob_pop[i]);
      fgets(aa,1000,fp);

      // if we're doing the FISHER mocks, scale mas by OMEGA_M and velocity by OMEGA_Z and GAMMA
      if(MCMC==3)
	{
	  mass_glob_pop[i] *= (OMEGA_TEMP/0.30711);
	  vxhs[i] *= vfac;
	  vyhs[i] *= vfac;
	  vzhs[i] *= vfac;
	}

      fscanf(fdelta,"%f", &delta_glob_pop[i]);
      fgets(aa,1000,fdelta);
      //fmuh(mass_glob_pop[i]);
    }
  fclose(fp);
  fclose(fdelta);

  for (i = 0; i < nhalos_glob_pop; ++i)
    {
      mass = mass_glob_pop[i];
      //M_min_delta(mass, delta_glob_pop[i], 1);
      if(mass > HOD.M_max) continue;
      if(mass < HOD.M_low) continue;

      ncen=N_cen(mass);
      box_n_glob_pop += ncen;
      nsat = N_sat(mass);
      box_n_glob_pop += nsat;
    }

  box_n_glob_pop /= (BOX_SIZE*BOX_SIZE*BOX_SIZE);


  //what's the density in normal HOD?
  printf("DENSITY: %e\n",box_n_glob_pop);
  //exit(0);
  

  // normalise this delta dependence (only if we need to)
  if(ROOTNORM<-1)
    rootnorm = zbrent(find_norm, -HOD.f_env, HOD.f_env, 0.001);
  else
    rootnorm = ROOTNORM;
  ROOTNORM=rootnorm;
  fprintf(stdout,"poplate_simulation_hod_density2> normalization: %f\n", rootnorm);
  
  // now set the random seed
  fprintf(stdout,"poplate_simulation_hod_density2> iseed: %d\n", IDUM_MCMC);
  srand48(IDUM_MCMC);


  for(ii = 0; ii < nhalos; ++ii)
    {
      haloid++; //in case it's not iterated in the file

      mass = mass_glob_pop[ii];
      xh[0] = xhs[ii];
      xh[1] = yhs[ii];
      xh[2] = zhs[ii];
      vh[0] = vxhs[ii];
      vh[1] = vyhs[ii];
      vh[2] = vzhs[ii];

      if(mass > HOD.M_max) continue;

      if(mass < mlo) continue;

      for (i = 0; i < 3; ++i)
        {

          if(xh[i]<0)xh[i]+=BOX_SIZE;
          if(xh[i]>BOX_SIZE)xh[i]-=BOX_SIZE;
        }

      i1 = (int)(log10(mass)/0.1);
      halocnt[i1]++;	  

      HOD.delta = delta_glob_pop[ii];
      M_min_delta(mass, delta_glob_pop[ii], rootnorm);

      ncen=N_cen(mass);
      ngals += ncen;

      if(drand48()>ncen)goto SATELLITES;

      if(VBIAS_C>0)
	{
	  NFW_central_velocity(mass,vg,mag);
	  for(i=0;i<3;++i)
	    vh[i]+=vg[i];
	}
      fprintf(fp2,"%.3f %.3f %.3f %.1f %.1f %.1f %d %.2f %.2f\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],0, log10(mass), delta_glob_pop[ii]);
      //fprintf(fpt,"%d\n",haloid);
 
      // THis mocks up the ~22 columns of the MR mocks
      //      fprintf(fpt,"%d %d %d %e %e %e 0.0 0.0 0.0 %e %e %e 0.0 0.0 0.0 0.0 0.0 0.0 %e 0.0 0.0\n",
      //	      0,0,0,xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],log10(mass));

      if(VBIAS_C>0)
	{
	  for(i=0;i<3;++i)
	    vh[i]-=vg[i];
	}
      

      //muh(i);
      //fmuh(ncen);
      //fmuh(mass);
      galcnt[i1]++;

      if(TRACK_GALAXIES)
	{
	  id2++;
	  galid[id2] = haloid;
	  galarr[id2][0] = xh[0];
	  galarr[id2][1] = xh[1];
	  galarr[id2][2] = xh[2];
	  galarr[id2][3] = vh[0];
	  galarr[id2][4] = vh[1];
	  galarr[id2][5] = vh[2];
	}

      if(TRACK_GALAXIES)
	{
	  id1++;
	  haloarr[id1][0] = xh[0];
	  haloarr[id1][1] = xh[1];
	  haloarr[id1][2] = xh[2];
	  haloarr[id1][3] = vh[0];
	  haloarr[id1][4] = vh[1];
	  haloarr[id1][5] = vh[2];
	  haloarr[id1][6] = mass;
	}

    SATELLITES:
      nsat = N_sat(mass);
      ngals += nsat;

      if(nsat>250)
	n1 = gasdev(&IDUM3)*sqrt(nsat) + nsat;
      else
	n1 = poisson_deviate(nsat);      
      /*
      if(nsat>1)
	fprintf(stdout,"BUH %d %f %e %d %e %f\n",haloid,nsat,mass,n1,HOD.M1,pow(mass/HOD.M1,HOD.alpha));
      */
      
      for(i=1;i<=n1;++i)
	{
	  r = NFW_position(mass,xg);
	  if(JEANS_DISPERSION)
	    jeans_dispersion(mass,r,vg);
	  else
	    NFW_velocity(mass,vg,mag);
	  for(k=0;k<3;++k)
	    {
	      xg[k]+=xh[k];
	      if(xg[k]<0)xg[k]+=BOX_SIZE;
	      if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
	      vg[k]+=vh[k];
	      /* vg[k] = vgf[k]; */
	      /*
		vg[k] = gasdev(&IDUM3)*500;
		vg[k] = non_gaussian_velocity();
	      */
	    }	
	  fprintf(fp2,"%.3f %.3f %.3f %.1f %.1f %.1f %d %.2f %.2f\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2],1, log10(mass), delta_glob_pop[ii]);
	  //fprintf(fpt,"%d\n",haloid);

	  //	  fprintf(fpt,"%d %d %d %e %e %e 0.0 0.0 0.0 %e %e %e 0.0 0.0 0.0 0.0 0.0 0.0 %e 0.0 0.0\n",1,1,1,xg[0],xg[1],xg[2],vg[0],vg[1],vg[2],log10(mass));

	  if(TRACK_GALAXIES)
	    {
	      id2++;
	      galid[id2] = haloid;
	      galarr[id2][0] = xg[0];
	      galarr[id2][1] = xg[1];
	      galarr[id2][2] = xg[2];
	      galarr[id2][3] = vg[0];
	      galarr[id2][4] = vg[1];
	      galarr[id2][5] = vg[2];
	    }

	  /* Bin up the galaxies by halo mass to check the HOD
	   */
	  galcnt[i1]++;	  
	}
    }
  fclose(fp2);
  //fclose(fpt);

  fprintf(stdout,"GALAXY_DENSITY:\t%e\nREAL DENSITY:\t%e\n",GALAXY_DENSITY, ngals/(BOX_SIZE*BOX_SIZE*BOX_SIZE));

  // reset the HOD to original values
  HOD.M_min = HOD.M_min_0;
  HOD.M_cut = HOD.M_cut_0;
  HOD.M_low = HOD.M_low_0;
  return;


  /* end ehere */

  /* output the binned HOD
   */
  sprintf(aa,"%s.binned_HOD",Task.root_filename);      
  fp2=fopen(aa,"w");
  for(i=0;i<1000;++i)
    if(galcnt[i]>0)
      fprintf(fp2,"%d %f %f %d %d\n",
	      i,(i+0.5)*0.1,(float)galcnt[i]/halocnt[i],galcnt[i],halocnt[i]);
  fclose(fp2);


  /* Calculate the two-halo term
   */
  if(TRACK_GALAXIES)
    {
      fprintf(stderr,"Calling nbody_2halo...\n");
      calc_nbody_two_halo(galarr,galid,id2);
    }
  return ;
  
  /* Track the close pairs - for Non-Tinkers, don't worry about this.
   */
  rmin = 1.5;
  rmax = 2.5;
  for(i=1;i<=id2;++i)
    for(j=i+1;j<=id2;++j)
      {
	if(galid[i]==galid[j])continue;
	dx = galarr[i][0] - galarr[j][0];
	if(dx>BOX_SIZE/2)dx = BOX_SIZE-dx;
	if(dx>rmax)continue;
	dy = galarr[i][1] - galarr[j][1];
	if(dy>BOX_SIZE/2)dy = BOX_SIZE-dy;
	if(dy>rmax)continue;
	dz = galarr[i][2] - galarr[j][2];
	if(dz>BOX_SIZE/2)dz = BOX_SIZE-dz;
	if(dz>rmax)continue;
	dr = sqrt(dx*dx+dy*dy+dz*dz);
	if(dr<rmin || dr>rmax)continue;

	i1 = galid[i];
	j1 = galid[j];
	dx = haloarr[i1][0] - haloarr[j1][0];
	if(dx>BOX_SIZE/2)dx = BOX_SIZE-dx;
	dy = haloarr[i1][1] - haloarr[j1][1];
	if(dy>BOX_SIZE/2)dy = BOX_SIZE-dy;
	dz = haloarr[i1][2] - haloarr[j1][2];
	if(dz>BOX_SIZE/2)dz = BOX_SIZE-dz;
	drh = sqrt(dx*dx+dy*dy+dz*dz);
	
	rv1 = pow(3*haloarr[i1][6]/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
	rv2 = pow(3*haloarr[j1][6]/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
	if(dr < rv1+rv2)
	  {
	    dx = galarr[i][3] - galarr[j][3];
	    dy = galarr[i][4] - galarr[j][4];
	    dz = galarr[i][5] - galarr[j][5];
	    printf("%e %e %e %e %e %e %e\n",dr,rv1+rv2,haloarr[i1][6],haloarr[j1][6],dx,dy,dz);
	    fflush(stdout);
	  }
      }

  exit(0);
}

