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

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

/* Internal function
 */
void calc_nbody_two_halo(float **gal, int *id, int ngal);
void M_min_delta(double m, double delta, double norm);
double find_norm(double norm);

double box_n_glob_pop;
double slope_glob_pop=-0.05;
float *mass_glob_pop, *delta_glob_pop;
int nhalos_glob_pop;

void M_min_delta(double m, double delta, double norm)
{
  double factor;
  factor = norm + slope_glob_pop * delta;
  if (factor > 0.05)
    HOD.M_min = HOD.M_min_0 * factor;
  else
    HOD.M_min = HOD.M_min_0 * 0.05;
}

double find_norm(double norm)
{
  float mass, ncen, ngals, nsat;
  int i;

  for (i = 0; i < nhalos_glob_pop; ++i)
    {
      mass = mass_glob_pop[i];
      if(mass > HOD.M_max) continue;
      if(mass < HOD.M_low) continue;

      M_min_delta(mass, delta_glob_pop[i], norm);
      //printf("%e %e %e\n",HOD.M_min_0, HOD.M_min,delta_glob_pop[i]);
      ncen=N_cen(mass);
      ngals += ncen;
      nsat = N_sat(mass);
      ngals += nsat;
      if(isnan(ngals)){ fmuh(mass); exit(0); }
    }
  ngals /= (BOX_SIZE*BOX_SIZE*BOX_SIZE);

  printf("\nNORMALISATION FACTOR:\t%f\nNUMBER DENSITY OF GALAXIES:\t%f\ndifference:\t%f\n", norm, ngals, ngals - box_n_glob_pop);
  //exit(0);
  return (ngals - box_n_glob_pop);
  
}

void populate_simulation_hod_density()
{
  FILE *fp, *fdelta, *fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt, *fgal;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  int nhalos, haloid;
  float *xhs, *yhs, *zhs, *vxhs, *vyhs, *vzhs, rootnorm;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv, ngals;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax;
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
  fdelta = openfile("input/halo_env_density.dat");

  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");
  sprintf(aa,"%s.mock_halo",Task.root_filename);      
  fpt = fopen(aa,"w");

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(i=0;i<1000;++i)
    galcnt[i]=0;

  set_HOD_params();
  HOD.M_min_0 = HOD.M_min;
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

  for(i = 0; i < nhalos; ++i)
    {
      fscanf(fp,"%f %f %f %f %f %f %f %f %f %f",
             &x1, &mass_glob_pop[i], &x1, &x1, &xhs[i], &yhs[i], &zhs[i], &vxhs[i], &vyhs[i], &vzhs[i]);
      fgets(aa,1000,fp);

      fscanf(fdelta,"%f", &delta_glob_pop[i]);
      fgets(aa,1000,fdelta);
    }
  fclose(fp);
  fclose(fdelta);

  for (i = 0; i < nhalos_glob_pop; ++i)
    {
      mass = mass_glob_pop[i];
      if(mass > HOD.M_max) continue;
      if(mass < HOD.M_low) continue;

      ncen=N_cen(mass);
      box_n_glob_pop += ncen;
      nsat = N_sat(mass);
      box_n_glob_pop += nsat;
    }

  box_n_glob_pop /= (BOX_SIZE*BOX_SIZE*BOX_SIZE);

  find_norm(-70.0);
  find_norm(-60.0);
  find_norm(-50.0);
  find_norm(-40.0);
  find_norm(-30.0);
  find_norm(-25.0);
  find_norm(-20.0);
  find_norm(-15.0);
  find_norm(-10.0);
  find_norm(-5.0);
  find_norm(0.0);
  find_norm(5.0);
  find_norm(10.0);
  find_norm(15.0);
  find_norm(20.0);


  exit(0);

  // normalise this delta dependence
  rootnorm = zbrent(find_norm, 0.1, 25.0, 0.0001);
  fp = fopen("output/delta_func.params", "w");
  fprintf(fp, "SLOPE:\t%f\nINTERCEPT:\t%f\n", slope_glob_pop, rootnorm);
  fclose(fp);

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
//      HOD.M_min = M_min_delta(mass, delta[i]);
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
      fprintf(fp2,"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
      fprintf(fpt,"%d\n",haloid);
 
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
	  fprintf(fp2,"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	  fprintf(fpt,"%d\n",haloid);

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
  fclose(fpt);

  fprintf(stdout,"GALAXY_DENSITY:\t%e\nREAL DENSITY:\t%e\n",GALAXY_DENSITY, ngals/(BOX_SIZE*BOX_SIZE*BOX_SIZE));

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

  exit(0);

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

