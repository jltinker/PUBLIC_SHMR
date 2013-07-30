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
double generate_satellite_galaxy_tab(double mass, double mstarlow, int reset);
void boxwrap_galaxy(float xh[], double xg[]);

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
double ms_to_mhalo_inversion(double mass);
double distance_redshift(double z);
double red_central_fraction(double mass, double *a);



void populate_simulation()
{
  FILE *fp,*fpa[9],*fp2,*fp3,*fpb[9],*fpc[9],*fps[9],*fpt,*fpblue, *fpsub,*fp2r,*fp3r;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000];
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1,*ibluecen,nhalo;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax, fred,p;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  //static int nhalo, *used_halo, iuse_flag = 1;
  int ihalo=0, ii;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;
  
  // new for mcmc_lensing;
  double m1, m2, mfslope, rad, mstar, rproj, xs[20], vhalo, imag, minimum_halo_mass, 
    minimum_mstar, minimum_lgmstar;
  long IDUM=-555;
  float xx[20], xd[3], xv[3], *subi, *submass, *temp;


  fprintf(stderr,"\n\nPOPULATING SIMULATION WITH SHMR.\n");
  fprintf(stderr,    "--------------------------------\n\n");

  fp=openfile(Files.HaloFile);
  fprintf(stderr,"opening file: [%s]\n",Files.HaloFile);
  nhalo = filesize(fp);
  fprintf(stderr,"number of lines: [%d]\n",nhalo);

  ibluecen = ivector(1,nhalo);

  if(Files.UseStellarMassBinsClustering)
    input_stellar_mass_bins();
  minimum_mstar = pow(10.0,wpl.mstar_wplo[0][1]);
  minimum_lgmstar = wpl.mstar_wplo[0][1];
  fmuh(minimum_lgmstar);

  if(!COLOR)
    {
      sprintf(aa,"%s.shmr_mock",Task.root_filename);      
      fp2 = fopen(aa,"w");
      fprintf(stderr,"outfile: [%s]\n",aa);
      sprintf(aa,"%s.shmr_mock_haloid",Task.root_filename);      
      fp3 = fopen(aa,"w");
    }
  else
    {
      sprintf(aa,"%s.shmr_mock_blue",Task.root_filename);      
      fp2 = fopen(aa,"w");
      fprintf(stderr,"outfile: [%s]\n",aa);
      sprintf(aa,"%s.shmr_mock_haloid_blue",Task.root_filename);      
      fp3 = fopen(aa,"w");
    }

  BLUE_FLAG = 1;
  wpl.reset_inversion = 1;
  set_up_hod_for_shmr(minimum_mstar,5.0E12,wpl.a);
  minimum_halo_mass = HOD.M_low; 
  generate_satellite_galaxy_tab(1.0E12,minimum_lgmstar,1);
  fprintf(stderr,"Expected number density of galaxies N(>%.2f Msol)= %e\n",log10(minimum_mstar),GALAXY_DENSITY);
  fprintf(stderr,"Minimum halo mass (M_low) required= %e\n",minimum_halo_mass);

  haloid = 0;
  for(ii=1;ii<=nhalo;++ii)
    {
      ibluecen[ii] = 0;
      fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
	     &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);

      // output the central with probability given by red_central_fraction
      fred = red_central_fraction(mass,wpl.a);
      p = drand48();
      if(p>fred)ibluecen[ii] = 1; // mark as blue or red before cutting for min halo mass
      if(mass<minimum_halo_mass)continue;
      haloid++;

      // decide whether this halo's center is red or blue
      if(p>fred) 
	{
	  // get the mass of the central galaxy
	  mstar = log10(ms_to_mhalo_inversion(mass));
	  m1 = gasdev(&IDUM)*wpl.a[6+(1-BLUE_FLAG)*11]+(mstar);
	  if(m1>minimum_lgmstar) {
	    fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],1,m1,mass);
	    fprintf(fp3,"%d\n",haloid);
	  }
	}

      // tabulate the CLF for this halo mass
      //nsat = generate_satellite_galaxy(-mass,9.0);
      nsat = generate_satellite_galaxy_tab(-mass,minimum_lgmstar,0);
      //printf("nsat check> %e\n",nsat);
      n1 = poisson_deviate(nsat);

      for(i=1;i<=n1;++i)
	{
	  //m1 = generate_satellite_galaxy(mass,9.0);
	  m1 = generate_satellite_galaxy_tab(mass,minimum_lgmstar,0);
	  r = NFW_position(mass,xg);
	  NFW_velocity(mass,vg,mag);
	  boxwrap_galaxy(xh,xg);
	  fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",
		  xg[0],xg[1],xg[2],vh[0]+vg[0],vh[1]+vg[1],vh[2]+vg[2],0,m1,mass);
	  fprintf(fp3,"%d\n",haloid);
	}	
      fflush(fp2);
      fflush(fp3);
    }
  fclose(fp2);
  fclose(fp3);

  /* ---------------------------------------------
   *  now got back and do the REDS, if required
   * ---------------------------------------------
   */
  if(!COLOR)return;

  sprintf(aa,"%s.shmr_mock_red",Task.root_filename);      
  fp2 = fopen(aa,"w");
  fprintf(stderr,"outfile: [%s]\n",aa);
  sprintf(aa,"%s.shmr_mock_haloid_red",Task.root_filename);      
  fp3 = fopen(aa,"w");
  
  rewind(fp);
  BLUE_FLAG = 0;
  wpl.reset_inversion = 1;
  set_up_hod_for_shmr(2.0E9,5.0E12,wpl.a);
  set_up_hod_for_shmr(1.0E9,5.0E12,wpl.a);
  minimum_halo_mass = HOD.M_low;
  fprintf(stderr,"Expected number density of red galaxies N(>10^9 Msol)= %e\n",GALAXY_DENSITY);
  fprintf(stderr,"Minimum halo mass (M_low) required= %e\n",minimum_halo_mass);

  generate_satellite_galaxy_tab(1.0E12,minimum_lgmstar,1);
  haloid = 0;
  for(ii=1;ii<=nhalo;++ii)
    {
      fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
	     &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);
      if(mass<minimum_halo_mass)continue;
      haloid++;

      // we already know whether this halo is blue central
      if(ibluecen[ii])goto POP_SATS;

      // get the mass of the central galaxy
      mstar = log10(ms_to_mhalo_inversion(mass));
      m1 = gasdev(&IDUM)*wpl.a[6+(1-BLUE_FLAG)*11]+(mstar);
      if(m1>minimum_lgmstar) {
	fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],1,m1,mass);
	fprintf(fp3,"%d\n",haloid);
      }

    POP_SATS:
      // tabulate the CLF for this halo mass
      nsat = generate_satellite_galaxy_tab(-mass,minimum_lgmstar,0);
      n1 = poisson_deviate(nsat);

      for(i=1;i<=n1;++i)
	{
	  m1 = generate_satellite_galaxy_tab(mass,minimum_lgmstar,0);
	  r = NFW_position(mass,xg);
	  NFW_velocity(mass,vg,mag);	
	  boxwrap_galaxy(xh,xg);
	  fprintf(fp2,"%e %e %e %e %e %e %d %e %e\n",
		  xg[0],xg[1],xg[2],vh[0]+vg[0],vh[1]+vg[1],vh[2]+vg[2],0,m1,mass);
	  fprintf(fp3,"%d\n",haloid);
	}	
      fflush(fp2);
      fflush(fp3);
    }

  fclose(fp);

  return ;
  
}

void boxwrap_galaxy(float xh[], double xg[])
{
  int i;
  for(i=0;i<3;++i)
    {
      xg[i] += xh[i];
      if(xg[i]>BOX_SIZE)xg[i]-=BOX_SIZE;
      if(xg[i]<0)xg[i]+=BOX_SIZE;
    }
}

/* Generate a random integer based on a Poisson distribution 
 * with mean given as input.
 */
int poisson_deviate(double nave)
{
  static int flag=0;
  double p,pp;
  int n;

  p=0;
  pp=1;

  while(p<pp)
    {
      if(nave<1)
	n=(int)(drand48()*20);
      else
	n=(int)(drand48()*30*nave);
      p=poisson_prob(n,nave);
      pp=drand48();
    }
  return(n);
}

/* Poisson probability of n given n_average
 */
double poisson_prob(int n, double nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((float)(fac*exp(-nave)));
}

/* Randomy generates a position away from the origin with 
 * a probability given by the NFW profile for a halo of the input
 * mass (and including the CVIR_FAC)
 */
double NFW_position(double mass, double x[])
{
  double r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir;
  
  cvir=halo_concentration(mass)*CVIR_FAC;
  rvir=pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=drand48()*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
    
    if(drand48()<=pr)
      {
	costheta=2.*(drand48()-.5);
	sintheta=sqrt(1.-costheta*costheta);
	signs=2.*(drand48()-.5);
	costheta=signs*costheta/fabs(signs);
	phi1=2.0*PI*drand48();
	
	x[0]=r*sintheta*cos(phi1);
	x[1]=r*sintheta*sin(phi1);
	x[2]=r*costheta;
	return r;
      }
  }
}

/* This is the NFW density profile
 */
double NFW_density(double r, double rs, double ps)
{
  return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
    fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19*sqrt(1+REDSHIFT);
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS;
  return(0);
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_central_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS_C;
  return(0);
}


double generate_satellite_galaxy(double mass, double mstarlow)
{
  static double *mgal, *ngal, *zz, prev_mass = -1;
  static int n=0;

  int i, j, k, iter=0;
  double mlo, mhi, dm, ngaltot, pmax, p, x, m;

  if(mass != prev_mass)
    {
      if(mass<0)mass = -mass;
      prev_mass = mass;
      if(!n)
	{
	  n = 20;
	  mgal = dvector(1,n);
	  ngal = dvector(1,n);
	  zz = dvector(1,n);
	}
      
      mlo = mstarlow;
      mhi = 12.0;
      dm = (mhi-mlo)/(n-1);

      ngaltot = 0;
      pmax = 0;
      for(i=1;i<=n;++i)
	{
	  mgal[i] = (i-1)*dm + mlo;
	  set_up_hod_for_shmr(pow(10.0,mgal[i]-dm/2),pow(10.0,mgal[i]+dm/2),wpl.a);
	  ngal[i] = N_sat(mass);
	  //printf("%e %f %f\n",mass,mgal[i],ngal[i]);
	  if(ngal[i]>pmax)pmax=ngal[i];
	}
      set_up_hod_for_shmr(pow(10.0,mstarlow),pow(10.0,12.0),wpl.a);    
      ngaltot = N_sat(mass);
      for(i=1;i<=n;++i)
	ngal[i] /= pmax;
      spline(mgal,ngal,n,1.0E+30,1.0E+30,zz);

      // if reset call, send total number of satellites
      return ngaltot;
    }

  // test
  while(iter<1000)
    {
      iter++;
      m = drand48()*(12.0-mstarlow)+mstarlow;
      splint(mgal,ngal,zz,n,m,&p);
      printf("TEST2 %f %e\n",m,p);
    }
  exit(0);

  // randomly grab
  while(1)
    {
      m = drand48()*(12.0-mstarlow)+mstarlow;
      splint(mgal,ngal,zz,n,m,&p);
      if(drand48()<p)break;
    }
  return m;
}


double generate_satellite_galaxy_tab(double mass, double mstarlow, int reset)
{
  static double **mgal, **ngal, **zz, prev_mass = -1, *mhalo, *totsat, *zzh, dlogmh,
    masslimit;
  static int n=0, nm, *nn;

  int i, j, k, iter=0;
  double mlo, mhi, dm, ngaltot, pmax, p, x, m, mhlo, mhhi, p1, w;

  if(reset)
    {
      if(!n)
	{
	  n = 100;
	  nm = 100;
	  mgal = dmatrix(1,nm,1,n);
	  ngal = dmatrix(1,nm,1,n);
	  zz = dmatrix(1,nm,1,n);
	  mhalo = dvector(1,nm);
	  totsat = dvector(1,nm);
	  zzh = dvector(1,nm);
	  nn = ivector(1,nm);
	}
      
      mlo = mstarlow;
      mhi = 12.0;
      dm = (mhi-mlo)/(n-1);
      mhlo = 1.0E+11;
      mhhi = 1.0E+16;
      dlogmh = log(mhhi/mhlo)/(nm-1);

      masslimit = mhlo;
      for(j=1;j<=nm;++j)
	{
	  ngaltot = 0;
	  pmax = 0;
	  mhalo[j] = exp(dlogmh*(j-1))*mhlo;
	  nn[j] = n;
	  for(i=1;i<=n;++i)
	    {
	      mgal[j][i] = (i-1)*dm + mlo;
	      set_up_hod_for_shmr(pow(10.0,mgal[j][i]-dm/2),pow(10.0,mgal[j][i]+dm/2),wpl.a);
	      p = N_sat(mhalo[j]);
	      if(i==1 && p<=0)masslimit = mhalo[j];
	      //if(BLUE_FLAG==0)
	      //printf("BUH %d %d %e %e %e\n",j,i,mhalo[j],mgal[j][i],p);
	      if(p<1.0E-10) { nn[j] = i-1; break; }
	      ngal[j][i] = p;
	      if(ngal[j][i]>pmax)pmax=ngal[j][i];
	    }
	  
	  set_up_hod_for_shmr(pow(10.0,mstarlow),pow(10.0,12.0),wpl.a);    	  
	  totsat[j] = log(N_sat(mhalo[j]));
	  //if(BLUE_FLAG==0)
	  //printf("BUHTOT %e\n",N_sat(mhalo[j]));
	  mhalo[j] = log(mhalo[j]);
	  //printf("> %e %e\n",mhalo[j]/log(10),totsat[j]/log(10));
	  for(i=1;i<=nn[j];++i)
	    ngal[j][i] /= pmax;
	  for(i=1;i<=nn[j];++i)
	    ngal[j][i] = log(ngal[j][i]);
	  spline(mgal[j],ngal[j],nn[j],1.0E+30,1.0E+30,zz[j]);
	}
      spline(mhalo, totsat,nm,1.0E+30,1.0+30,zzh);
      return;
    }

  //if mass<0 just return total number of satellites
  if(mass<0)
    {
      if(-mass<masslimit)return 0;
      splint(mhalo, totsat, zzh, nm,log(-mass),&p);
      return exp(p);
    }

  // find the mass bins the bracket the input mass
  mass = log(mass);
  for(j=2;j<=nm;++j)
    if(mhalo[j]>mass)break;
  if(j>nm)j=nm;
  w = (mass-mhalo[j-1])/dlogmh;

  // randomly grab 
  while(1)
    {
      iter++;
      m = drand48()*(12.0-mstarlow)+mstarlow;
      splint(mgal[j],ngal[j],zz[j],nn[j],m,&p);
      p = exp(p);
      if(m>mgal[j][nn[j]])p = 0;
      splint(mgal[j-1],ngal[j-1],zz[j-1],nn[j-1],m,&p1);
      p1 = exp(p1);
      if(m>mgal[j-1][nn[j-1]])p1 = 0;
      p = p*w + p1*(1-w);
      if(iter>10000){
	printf("BOO %e %e %e %e %e\n",exp(mass),m,p,w,p1); m = 9.0; break; }
      if(drand48()<p)break;
    }
  return m;
}


