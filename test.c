#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

void deconvolved_smf(void);
void bias_from_simulation(void);
int optimal_spacing();
void latin_hypercube(int);
double distance_redshift(double);
double boss_redshift_evolution(void);

void test(int argc, char **argv)
{
  int i,i2,i1,n,j;
  float xx[10],x1,x2,x3,x4,x5,r,dr,x0;
  double xp[3],vv[3];
  double m, x, density, r1, r2, volume, bb, a;
  FILE *fp, *fp1;
  char aa[100], fname[100];

  // print out random Poisson deviates for every halo in the box
  fp = openfile("nsat.dat");
  n = filesize(fp);
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%f %f",&x1,&x2);
      printf("%e %e %d\n",x1,x2,poisson_deviate(x2));
    }
  exit(0);

  shambias_loop();
  exit(0);



  // print out a bunch of random NFW positions and velocities
  for(i=1;i<=10000;++i)
    {
      NFW_position(1.0E12,xp);
      NFW_velocity(1.0E+12,vv,-1);
      printf("NFW %f %f %f %f %f %f\n",xp[0],xp[1],xp[2],vv[0],vv[1],vv[2]);
    }
  exit(0);


  // get the spergel issue: Mavg=1.9E14 Msol/h bias ratio = 2

  // standard HOD
  HOD.M_min = 0;
  GALAXY_DENSITY = 0.003;
  HOD.M1 = pow(10.0,13.46);
  HOD.M_cut = pow(10.0,12.60);
  HOD.alpha = 1.03;
  HOD.sigma_logM = 0.3;

  set_HOD_params();
  x1 = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  x2 = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f %f\n",x1,x2);

  // normalizing
  ASSEMBLY_BIAS = 0;
  x1 = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = -1;
  x2 = 1./2.*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = 1;
  x3 = 1./2.*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f\n",x1+x2+x3);
  x0 = x1 + x2 + x3;


  //maximal issue: what if all galaxies are in the high-biased halos.
  ASSEMBLY_BIAS = 1;
  x1 = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  x2 = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f\n",x1/x0-1);

  // what if all the satellites are in the high-biased halos
  ASSEMBLY_BIAS = 0;
  x1 = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = 1;
  x2 = qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f\n",(x1+x2)/x0-1);


  // what if twice as many satellites in the high-biased halos
  ASSEMBLY_BIAS = 0;
  x1 = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = -1;
  x2 = 1./3.*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = 1;
  x3 = 2./3.*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f\n",(x1+x2+x3)/x0-1);

  // what if twice as many satellites in the high-biased halos
  ASSEMBLY_BIAS = 0;
  x1 = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = -1;
  x2 = 1./2.2*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  ASSEMBLY_BIAS = 1;
  x3 = 1.2/2.2*qromo(func_satellite_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"%f\n",(x1+x2+x3)/x0-1);


  exit(0);

  // keep nbar the same. keep fsat the same.
  


  // output matter correlation function and power spectrum as a function of expansion parameter
  OMEGA_M = 0.3;
  OMEGA_L = 0.7;
  for(i=1;i<=10;++i)
    {
      a = i/10.;
      REDSHIFT = 1/a - 1;
      SIGMA_8 = 0.8*growthfactor(REDSHIFT);
      RESET_COSMOLOGY++;
      sprintf(Task.root_filename,"lcdm%d",i);
      output_matter_correlation_function();
      //output_matter_power_spectrum();
    }
  exit(0);



  //  HIGH_PRECISION = 1;
  for(i=1;i<=10;++i)
    {
      a = i/10.;
      REDSHIFT = 1/a - 1;
      SIGMA_8 = 0.8*growthfactor(REDSHIFT);
      RESET_COSMOLOGY++;
      sprintf(Task.root_filename,"mfi%d",i);
      output_halo_mass_function();
      output_matter_variance();
    }
  exit(0);


  // inputs from the cosmomc chains
  fp = openfile("pk_inputs_ow0wacdm_planck_bao_snu.dat");
  n = filesize(fp);
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&SPECTRAL_INDX, &OMEGA_L, &OMEGA_M, &SIGMA_8, &HUBBLE, &OMEGA_B);
      RESET_COSMOLOGY++;
      sprintf(Task.root_filename, "pk_%d", i);
      output_matter_power_spectrum();
    }
  exit(0);

  for(i=1;i<=100;++i)
    {
      OMEGA_M = 0.3;
      OMEGA_L = 0.7;
      x1 = growthfactor(1/(i/100.)-1);
      OMEGA_L = 0.0;
      x2 = growthfactor(1/(i/100.)-1);
      printf("%f %f %f\n",1/(i/100.)-1,x1,x2);
    }
  exit(0);
      

  boss_redshift_evolution();
  exit(0);

  //get volume between z=1 and 2
  r1 = distance_redshift(1.0);
  r2 = distance_redshift(2.0);
  fmuh(r1);
  fmuh(r2);
  volume = 4./3.*PI*(r2*r2*r2 - r1*r1*r1)*7500./41253.;
  density = 15*7500./volume;
  fmuh(density);

  // okay, let's find an HOD model for this.
  HOD.pdfs = 0;
  HOD.M_min = 0;
  HOD.MaxCen = 0.01;
  HOD.sigma_logM = 0.5;
  GALAXY_DENSITY = density;
  set_HOD_params();
  bb = qromo(func_galaxy_bias, log(HOD.M_low), log(HOD.M_max), midpnt)/GALAXY_DENSITY;

  printf("BIAS: %f %e %e\n",bb, HOD.M_min, HOD.M_low);

  exit(0);

  for(i=0;i<=100;++i)
    printf("%e %e\n",i/100.,growthfactor(i/100.0));
  exit(0);

  latin_hypercube(7);
  // latin_hypercube(8);
  exit(0);
  latin_hypercube(2);
  latin_hypercube(1);
  exit(0);
  latin_hypercube(2);
  exit(0);
  latin_hypercube(5);
  exit(0);
  latin_hypercube(3);
  exit(0);
  latin_hypercube(4);
  exit(0);

  latin_hypercube(-2);
  exit(0);

  optimal_spacing();
  exit(0);


  fp = openfile(argv[3]);
  fscanf(fp,"%lf %lf %lf %lf %lf %lf",&x,&HOD.M_min,&HOD.M1,&HOD.alpha,&HOD.M_cut,&HOD.sigma_logM);
  GALAXY_DENSITY = 0;
  set_HOD_params();
  fmuh(GALAXY_DENSITY);

  RESOLUTION = 1.0;
  m = OMEGA_M*RHO_CRIT*100*pow(RESOLUTION,3.0);
  x = qromo(func_galaxy_density,log(m),log(HOD.M_max),midpnt);
  fprintf(stderr,"res: %f min_halo_mass: %e frac_in_halos: %f\n",RESOLUTION,m,x/GALAXY_DENSITY);

  RESOLUTION = 0.8;
  m = OMEGA_M*RHO_CRIT*100*pow(RESOLUTION,3.0);
  x = qromo(func_galaxy_density,log(m),log(HOD.M_max),midpnt);
  fprintf(stderr,"res: %f min_halo_mass: %e frac_in_halos: %f\n",RESOLUTION,m,x/GALAXY_DENSITY);

  if(argc>4)
    {
      RESOLUTION = atof(argv[4]);
      m = OMEGA_M*RHO_CRIT*100*pow(RESOLUTION,3.0);
      x = qromo(func_galaxy_density,log(m),log(HOD.M_max),midpnt);
      fprintf(stderr,"res: %f min_halo_mass: %e frac_in_halos: %f\n",RESOLUTION,m,x/GALAXY_DENSITY);
    }
  exit(0);


  populate_simulation_hod_density();

  fp = fopen("output/analytic_nbar", "w");
  fprintf(fp, "%e\n", GALAXY_DENSITY);
  fclose(fp);
  exit(0);

  //fisher();
  shambias_loop();
  bias_from_simulation();
  deconvolved_smf();

  // what's the fraction of galaxies resolved in halos of mass M 
  set_HOD_params();
  for(i=110;i<=160;++i)
    {
      m = pow(10.0,i/10.0);
      x = 1200.*pow(m/(100*OMEGA_M*RHO_CRIT),THIRD);
      printf("%d %e %f\n",i,qromo(func_galaxy_density,
			       i/10.0*log(10),log(HOD.M_max),midpnt)/GALAXY_DENSITY,x);
    }
  exit(0);

}
