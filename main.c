#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* test file routines.
 */
void test(int argc, char **argv);
double find_maximum_halo_mass(void);

int main(int argc, char **argv)
{
  double s1, delta_vir, omega_m, x;
  int i, j;
  FILE *fp;

#ifdef PARALLEL
  printf("STARTING>>>\n");
  fflush(stdout);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  printf("TASK %d reporting for duty.\n",ThisTask);
  fflush(stdout);
#endif

  ARGC = argc;
  ARGV = argv;

  OUTPUT=0;
  HOD.fredc = HOD.freds = 1.0;

  for(i=1;i<=99;++i)
    HOD.free[i]=0;
  wp.esys=0;

  Work.chi2=0;
  Work.imodel=1;

  USE_ERRORS = 0;
  ITRANS=4;
  HUBBLE=0.7;
  BEST_FIT = 0;
  HOD.M_sat_break = 1.0e14;
  HOD.alpha1 = 1.0;

  if(argc==1)
    endrun("./HOD.x hod.bat_file > output");

  read_parameter_file(argv[1]);

  if(REDSHIFT>0)
    {
      SIGMA_8 = SIGMA_8*growthfactor(REDSHIFT);
      HUBBLEZ = sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+1-OMEGA_M);
      OMEGA_Z = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
      fprintf(stdout,"SIGMA_8(Z=%.3f)= %.4f\n",REDSHIFT,SIGMA_8);
      fprintf(stdout,"H(Z=%.3f)/H0= %.4f\n",REDSHIFT,HUBBLEZ);
      HOD.M_min = 0;
      RESET_COSMOLOGY++;
      set_HOD_params();
    }

  /* Output the virial overdensity for reference.
   */
  if(OUTPUT)
    {
      omega_m=OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
      x=omega_m-1;
      delta_vir=(18*PI*PI+82*x-39*x*x)/(1+x);
      printf("DELTA_VIR(Omega_m,z) = %f\n",delta_vir);
    }

  /* Do some initialization if we're doing SHMR
   */
  if(SHMR_FLAG)
    {
      if(SATELLITE_PARAMETERIZATION)SHMR_PARAMS = 14;
      if(VARIABLE_ALPHA)SHMR_PARAMS += 2;
      if(VARIABLE_EXCLUSION)wpl.a[SHMR_PARAMS+1] = EXCLUSION_RADIUS;
      wpl.ncf = SHMR_PARAMS + VARIABLE_EXCLUSION;

      HOD.pdfs = 100;      
      HOD.pdfc = 101;
      wpx.calculate_two_halo = 1;
      input_stellar_mass_bins();
      // if we have input from the prompt, take that
      if(argc>2 && atoi(argv[2])!=999)
	{
	  fp = openfile(argv[2]);
	  fscanf(fp,"%d %d",&i,&j);
	  for(i=1;i<=wpl.ncf;++i)
	    fscanf(fp,"%lf",&wpl.a[i]);
	  fclose(fp);
	}	  
    }

  for(i=1;i<=wpl.ncf;++i)
    printf("wpl.a[%d]= %e\n",i,wpl.a[i]);

  /* LENSING TESTING FOR ALEXIE
   */
  if(argc>2)
    IDUM_MCMC=atoi(argv[2]);
  SIGMA_8Z0 = 0.8;

  if(argc>2)
    if(atoi(argv[2])==999)
      test(argc,argv);

  /* If there's no cross-correlation function,
   * set the second number density equal to the first
   */
  if(!XCORR)
    GALAXY_DENSITY2 = GALAXY_DENSITY;

  /* Initialize the non-linear power spectrum.
   */
  nonlinear_sigmac(8.0);
  sigmac_interp(1.0E13);
  sigmac_radius_interp(1.0);

  /* Skip the HOD stuff if we're SHMR-ing it:
   */
  if(SHMR_FLAG)
    {
      if(argc>2 && atoi(argv[2])==999)test(argc, argv);
      if(argc>3 && atoi(argv[3])==999)test(argc, argv);
      goto TASKS;
    }

  /* Get the galaxy bias factor
   */
  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"Galaxy Bias bg= %f\n",GALAXY_BIAS);
  fflush(stdout);

  /* Get the galaxy satellite fraction
   */
  s1=qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/
    GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"fsat %e\n",s1);
  fflush(stdout);

  /* Mean halo mass.
   */
  if(OUTPUT)
    fprintf(stdout,"M_eff %e\n",number_weighted_halo_mass());
  fflush(stdout);

  /* Set up BETA for wp integration.
   */
  BETA = pow(OMEGA_M,0.6)/GALAXY_BIAS;
  if(OUTPUT)
    printf("BETA = %f\n",BETA);

 TASKS:
  tasks(argc,argv);
}

/* This pair of functions is meant to find the halo mass at which 1 halo
 * should exist in a volume of size BOX_SIZE^3
 */
double func_max_mass(double mlo)
{
  return qromo(func_halo_density,mlo,log(1.0E17),midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE - 1;
}
double find_maximum_halo_mass()
{
  return zbrent(func_max_mass,log(1.0E12),log(1.0E16),1.0E-4);
}
