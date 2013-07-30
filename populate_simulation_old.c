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
double NFW_central_velocity(double mass, double v[], double mag);
void calc_nbody_two_halo(float **gal, int *id, int ngal);
void calc_nbody_one_halo(float **gal, int *id, int ngal);
double generate_satellite_galaxy(double mass, double mstarlow);
double random_magnitude(int ngal, float *mass, float *mag, double mgal);

double ms_to_mhalo_inversion(double mass);
double red_central_fraction(double mass, double *a);



void populate_simulation_old()
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt,*fpblue, *fpsub;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000];
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax, fred;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  //static int nhalo, *used_halo, iuse_flag = 1;
  int ihalo=0;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;
  
  // new for mcmc_lensing;
  double m1, m2, mfslope, rad, mstar, rproj, xs[20], vhalo, imag;
  long IDUM=-555;
  float xx[20], xd[3], xv[3], *subi, *submass, *temp;
  int *ii;

  // tabulte the satellite function fot Mgal=8.8
  double *mtest, *ntest, *ytest, mhi, dlogm, a;
  int nn=100;

  mtest=dvector(1,nn);
  ntest=dvector(1,nn);
  ytest=dvector(1,nn);

  mlo = 1.0E+10;
  mhi = 1.0E+15;
  dlogm = log(mhi/mlo)/(nn-1);
  set_up_hod_for_shmr(pow(10.0,8.8),pow(10.0,12.0),wpl.a);
  for(i=1;i<=nn;++i)
    {
      mtest[i] = exp((i-1)*dlogm)*mlo;
      ntest[i] = N_sat(mtest[i]);
    }
  spline(mtest,ntest,nn,1.0E+30,1.0E+30,ytest);


  fp=openfile(Files.HaloFile);
  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");

  fprintf(stderr,"opening file: [%s]\n",Files.HaloFile);
  fprintf(stderr,"outfile: [%s]\n",aa);


  while(!feof(fp))
    {
      // NB for L120b file:
      fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
	     &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);

      // correct for 30% error in fof_gadget mass
      //mass /= 1.15;

      // get the mass of the central galaxy
      mstar = log10(ms_to_mhalo_inversion(mass));
      m1 = gasdev(&IDUM)*wpl.a[6+(1-BLUE_FLAG)*11]+(mstar);
      
      // output the central with probability given by red_central_fraction
      fred = red_central_fraction(mass,wpl.a);
      if(BLUE_FLAG) fred = 1-fred;
      if(drand48()<=fred && m1>8.8)
	fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],1,m1,mass);

      splint(mtest,ntest,ytest,nn,mass,&nsat);
      if(nsat<1.0E-3)continue;

      // tabulate the CLF for this halo mass
      nsat = generate_satellite_galaxy(-mass,8.8);
      //nsat = generate_satellite_galaxy(-mass,7.0); // <-- for josh
      n1 = poisson_deviate(nsat);

      for(i=1;i<=n1;++i)
	{
	  m1 = generate_satellite_galaxy(mass,8.8);
	  r = NFW_position(mass,xg);
	  NFW_velocity(mass,vg,mag);
	  fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",
		  xh[0]+xg[0],xh[1]+xg[1],xh[2]+xg[2],vh[0]+vg[0],vh[1]+vg[1],vh[2]+vg[2],0,m1,mass);
	}	
      
      if(feof(fp))break;
    }
  fclose(fp2);
  fclose(fp);

  return ;
  
}
