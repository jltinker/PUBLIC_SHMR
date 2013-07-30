#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

double func_dr1(double z);
double func_wth1(double r);
double func_wth1_tabulated(double r);
double distance_redshift(double z);
double redshift_distance(double r);

struct WTH {
  double zmin, zmax;
  double rad, theta;
} wth;

double wtheta(double theta)
{
  double r1, r2, z1, z2, z;
  
  double rmin, rmax, dr, dz, Deltar, pz, s1, zmin, zmax, rad, drdz, omega_temp, s2;
  int nrad=200, i;

  
  static double *nquadri,*zquadri,*yquadri;
  double dz1,dz2;
  static int flag = 1, ntable;
  FILE *fp;

  // convert from arcsecs to radians!
  theta *= PI/180.0/3600.0;

  // implement with the tabulated N(z)
  if(flag)
    {
      flag = 0;
      // JPL
      //fp = openfile("/home/asleauth/HOD/Ccode/cosmos_redshift_distribution.dat");
      //fp = openfile("/Users/alexie/Work/HOD/Ccode/cosmos_redshift_distribution.dat");
      fp = openfile("cosmos_redshift_distribution.dat");
      
      ntable = filesize(fp);
      zquadri = dvector(1,ntable);
      nquadri = dvector(1,ntable);
      yquadri = dvector(1,ntable);
      
      for(i=1;i<=ntable;++i)
	{
	  fscanf(fp,"%lf %lf",&zquadri[i],&nquadri[i]);
	  nquadri[i] = log(nquadri[i]);
	}
      spline(zquadri,nquadri,ntable,1.0E+30,1.0E+30,yquadri);
    }
  
  /* for the COSMOS data
   */
  zmin = wpl.zlo;
  zmax = wpl.zhi;
  nrad = 20;

  wth.theta = theta;
  rmin = distance_redshift(zmin);
  rmax = distance_redshift(zmax);

  dr = (rmax - rmin)/nrad;
  Deltar = rmax - rmin;

  // actually, this is too far. truncate to 50 Mpc/h
  Deltar = 50.0;
  //Deltar = 30; // for risa (used in all calcs for paper)

  pz = 1.0/(zmax - zmin);
  pz *= pz;

  s1 = s2 = 0;
  

  //trying to implement charlie's way explicitly
  for(i=0;i<nrad;++i)
    {
      z = (zmax-zmin)/(nrad)*(i+0.5) + zmin;
      wth.rad = rad = distance_redshift(z);
      dz = (zmax-zmin)/nrad;
      //splint(zquadri,nquadri,yquadri,ntable,z,&pz);
      //pz = exp(pz)*exp(pz);
      pz = 1;
      drdz = c_on_H0/sqrt(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M));
      //drdz = c_on_H0*c_on_H0/(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M));
      s1 += 2*qromo(func_wth1,0.0,Deltar,midpnt)*pz/drdz*dz;
      s2 += sqrt(pz)*dz; //test
    }
  //printf("%e %e %e\n",s1,s2,dz);

  return s1/s2/s2;


  for(i=0;i<nrad;++i)
    {
      rad = (i-0.5)*dr + rmin;
      wth.rad = rad;
      z = redshift_distance(rad);
          
      splint(zquadri,nquadri,yquadri,ntable,z,&pz);
      pz = exp(pz)*exp(pz);

      /*
      SIGMA_8 = 0.9*growthfactor(z);
      RESET_COSMOLOGY++;
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      */
      //pz = exp(-(z-4.8)*(z-4.8)/(2*.4*.4));
      //pz = pz*pz;
      drdz = c_on_H0*c_on_H0/(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M));
      s1 += 2*qromo(func_wth1,0.0,Deltar,midpnt)*pz/drdz*dz;
      //s1 += 2*qtrap(func_wth1,0.0,Deltar,1.0E-5)*pz/drdz*dr;
      //s1 += 2*qromo(func_wth1_tabulated,0.0,Deltar,midpnt)*pz/drdz*dr;
      s2 += sqrt(pz/drdz)*dr;
		   
    }
  printf("%e %e %e\n",s1,s2,dz);

  return s1/s2/s2;
}

double func_wth1(double rlos)
{
  double r;
  r = sqrt(rlos*rlos + wth.theta*wth.theta*wth.rad*wth.rad);
  //return one_halo_real_space(r);
  //return two_halo_real_space(r);
  
  /*
  fprintf(stderr,"%e %e %e %e %e\n",HOD.M_min, HOD.M1, HOD.M_low, HOD.sigma_logM, GALAXY_DENSITY);
  fprintf(stderr,"1> %e\n",r);
  fprintf(stderr,"2> %e\n",one_halo_real_space(r));
  fprintf(stderr,"3> %e\n",two_halo_real_space(r));
  */
  return xi_interp(r);
  return one_halo_real_space(r) + two_halo_real_space(r);
  return one_halo_real_space(r) + linear_kaiser_distortion(r,rlos);
}

double func_wth1_tabulated(double rlos)
{
  float x1;
  char aa[1000];

  FILE *fp;
  double rad, xi;
  int i;
  static int flag = 1, n;
  static double *r, *x, *y;
  if(flag)
    {
      flag = 0;
      
      //fp = openfile("risa_xir.080.drg1.dat");
      //fp = openfile("xir.080.drg1x2.dat");
      fp = openfile("xitab.dat");
      n = filesize(fp);
      r = dvector(1,n);
      x = dvector(1,n);
      y = dvector(1,n);
     
      
      for(i=1;i<=n;++i)
	{
	  //fscanf(fp,"%f %f",&x1,&x1);
	  fscanf(fp,"%lf %lf %lf",&r[i],&x[i],&rad);
	  //fgets(aa,1000,fp);
	  if(x[i]<0){ n = i-1; break; }
	  r[i] = log(r[i]);
	  x[i] = log(x[i]);
	}
      spline(r,x,n,1.0E+30,1.0E+30,y);
    }

  rad = log(sqrt(rlos*rlos + wth.theta*wth.theta*wth.rad*wth.rad));
  splint(r,x,y,n,rad,&xi);
  return(exp(xi));
}

double redshift_distance(double r)
{
  int i;
  double z,dz,a;
  static double *rad,*redshift,*yy;
  static int flag=1, n=100;

  if(flag)
    {
      flag = 0;
      rad = dvector(1,n);
      redshift = dvector(1,n);
      yy = dvector(1,n);

      dz = 3.0/n;

      for(i=1;i<=n;++i)
	{
	  redshift[i] = 1 + dz*(i-1);
	  rad[i] = distance_redshift(redshift[i]);
	}
      spline(rad,redshift,n,1.0E+30,1.0E+30,yy);
    }
  splint(rad,redshift,yy,n,r,&a);
  return a;
}

double distance_redshift(double z)
{
  return c_on_H0*qtrap(func_dr1,0.0,z,1.0E-5);
}
double func_dr1(double z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}

void nden_check()
{
  double r1, r2, z1, z2, z, r, dz;
  
  double rmin, rmax, dr, Deltar, pz, x1, x2, zmin, zmax, rad, drdz, omega_temp;
  int nrad=20, i;

  zmin = 2.;
  zmax = 3;

  pz = 1.0/(zmax - zmin);
  pz *= pz;

  dz = (zmax-zmin)/nrad;

  x1 = x2 = 0;
  for(i=0;i<nrad;++i)
    {
      z = zmin + (i+0.5)*dz;
      SIGMA_8 = 0.9*growthfactor(z);
      RESET_COSMOLOGY++;
      rad = distance_redshift(z);
      drdz = c_on_H0/sqrt(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M))*rad*rad;
      //pz = exp(-(z-4)*(z-4)/(2*.3*.3));
      x1 += qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt)*drdz*pz;
      x2 += drdz*pz;
    }
  fmuh(x1/x2);
			     
  exit(0);
}
