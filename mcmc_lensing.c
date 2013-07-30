#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* External functions
 */
double ms_to_mhalo(double ms, double *a);  // in chi2_lensing.c, ms NOT in log
double generate_satellite_galaxy(double mass, double mstarlow); // from populate_simulation.c
void populate_simulation_old(void);

/* Internal functions.
 */
void wp_lensing_input(void);
double chi2_wp_lensing_wrapper(double *a);
double chi2_stellar_mass_function(double *a);
double mcmc_lensing_initialize(double *a, double **cov1, double *avg1, double *start_dev);
double func_satellite_density_xigm(double m);              
double funct_int_vary_scat(double ms);
double integral_constraint_correction(double r);
double red_central_fraction(double mass, double *a);
void fisher_matrix(double *a, double *eval, double **evect, int ishuffle);
void chain_post_process_single(double *a);
double chi2_red_central_fraction(double *a);
double chi2_lensing_total(double *a);
void initial_lensing_values(double *a, double **pp, double *yy);
void chi2_lensing_minimization();
void chain_post_process_clf(double *a);
int check_lensing_parameters(double *a);
void normalize_eigenvalues(double **evect, double *eval, double *a);
double func_eigen_normalize(double fac);
void chain_post_process(void);
void input_chain(double *eval, double **evect);
void chain_post_process_central_red_error(void);
double func_mean_hmass(double m);


/* Local global variable*/
double aa6,aa7,aa8,aa9,aa10,mh_scatter_vary;
double *eval_gl1, **evect_gl1, chi2original_gl1;
int eigeni_gl1;

/* FISHER CODES:
 * 2 - use fisher for covar
 * 1 - calc fisher matrix
 * 3 - renormalize previously calculated fisher matrix
 * 4 - use fisher to shuffle initial parmaeter set
 * 5 - input an external chain to get covariance matrix
 * 6 - chain post processing
 */

/* This is a flag to tell whether we're doing the 
 * red gal HOD or the blue gal HOD
 */
int MAKE_MOCKS = 0;
int ALPHA_FREE = 1;
int HOD_PARAMS = 10;
int CHI2_MINIMIZATION = 0;

/*----------------------------------------------------------------
 *        MCMC_LENSING
 * --------------------------------------------------------------*/

void mcmc_lensing()
{
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=3000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight,USE_IWEIGHT=0;
  double t0,tprev,temp,chi2a,chi2b,chi2wp,chi2lensing,chi2smf,r,chi2fred;
  double mh_test;
  int pcheck_count;
  double t1,t2,t3,t4;


  wpl.reset_fred = 1;
  wpl.mlow_flag = 0;
  MCMC=Task.MCMC;
  EXCLUSION = 4;    // Changes how 2h is calculated (4=mcmc_lensing, 5=chi2_lensing) 
  HOD.pdfc = 100;   // Functional form for Ncen (TRESHOLD SAMPLES)
  HOD.pdfs = 66;    // Functional form for Nsat (USED TO BE 6 -> CHANGED TO 66)

  //chain_post_process_single();


  pcheck=calloc(ptot,sizeof(int));
  pcheck_count = 0.0;

  /* SET THE COSMOLOGY
   * sigma8 is set in the hod.bat file,
   * as is REDSHIFT
   */
  SIGMA_8Z0 = SIGMA_8;
  SIGMA_8 = SIGMA_8Z0*growthfactor(REDSHIFT);
  RESET_COSMOLOGY++;

  // Change this to fit or not fit the lensing and SMF
  DONT_FIT_SMF        = 0;   // don't use the smf in the fit, set to 1 to not use SMF
  DONT_FIT_LENSING    = 0;   // don't use the lensing in the fit, set to 1 to not use lensing
  DONT_FIT_CLUSTERING = 0;   // set to 1 to not use the clustering

  DONT_FIT_SMF        = wpl.dont_fit_smf;   // don't use the smf in the fit, set to 1 to not use SMF
  DONT_FIT_LENSING    = wpl.dont_fit_lensing;   // don't use the lensing in the fit, set to 1 to not use lensing
  DONT_FIT_CLUSTERING = wpl.dont_fit_clustering;   // set to 1 to not use the clustering

  // read in the clustering data.
  if(MAKE_MOCKS==0)
    wp_lensing_input();

  /* SET RANDOM NUMBER GENERATOR
   * this doesn't set the parameter selection, which use gadev() 
   * but is used in drand48() calls
   */
  srand48(32498793);

  /* SET THE NUMBER OF FREE PARAMTERS
   */
  n = wpl.ncf = 10;   // With slopes for Mcut and Msat (default value)
  if(ALPHA_FREE) n = wpl.ncf = 11;   // With alpha, a[11]
  HOD_PARAMS = n; // this is the number of HOD params for EACH red/blue (without fredcen)

  /* QUENCHED: add new parameters: 
   * - 5 for the quenched central as function of mass (spline fit)
   * - 5 for new SHMR for red centrals (old parameters are for blue centrals)
   * - 4 for red satellites
   * - 1 for scatter in red centrals
   * -----
   *   15 new parameters, total of 25 parameters.
   */
  n = wpl.ncf = 2*n + 5;

  //if CHAIN processing, do it here
  if(FISHER==6)
    chain_post_process();
  //if CHAIN processing, do it here
  if(FISHER==7)
    chain_post_process_central_red_error();


  if(CHI2_MINIMIZATION)
    chi2_lensing_minimization();


  //** FOR PLOT: By hand **:
  // NB-- making this an input parameter
  // SET TO 2 FOR PRETTY PLOT
  //  LENSING_OUTPUT_FLAG = 0;

  if(OUTPUT)
    printf("mcmc_min> %d  free parameters\n",n);

  a=dvector(1,n);
  start_dev=dvector(1,n);
  aprev=dvector(1,n);
  atemp=dvector(1,n);
  cov1=dmatrix(1,n,1,n);
  avg1=dvector(1,n);

  tmp=dmatrix(1,n,1,n);
  tmp1=dmatrix(1,n,1,1);
  evect=dmatrix(1,n,1,n);
  eval=dvector(1,n);
  eval_prev=dvector(1,n);

  chain=dmatrix(1,imax_chain,1,n);
  iweight = ivector(1,imax_chain);
  for(i=1;i<=imax_chain;++i)
    iweight[i] = 0;

  IDUM=IDUM_MCMC;

  chi2prev=mcmc_lensing_initialize(a,cov1,avg1,start_dev);
  niter++;
  for(i=1;i<=n;++i)
    {
      aprev[i] = a[i];
      chain[1][i] = a[i];
    }

  if(FISHER && FISHER!=5){
    fisher_matrix(a,eval,evect,0);
    convergence=1;
    goto SKIP_BURN;
  }
  if(FISHER==5){
    input_chain(eval,evect);
    convergence=1;
    goto SKIP_BURN;
  }
  pcnt=0;
  pcheck[pcnt]=1;

  /* this is the initial loop at the beginning of the burn-in phase.
   * the first NSTEP elements in the chain are selected from hard-wired
   * dispersions, not using covariance matrix
   * This is the first NSTEP (~50) acceptances
   */
  //stepfac=1;
  stepfac = 0.2;

  while(niter<NSTEP)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}

      /* hardwire stepfac
       * Can tweak this here if the step size is getting too small at the
       * begining of the chain (or could also change STD_DEV)
       */
      stepfac=0.01;

      for(i=1;i<=n;++i)
	wpl.a[i] = a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];
      wpl.reset_inversion = 1; // reset the inversion of the ms-mh inversion function

      // Check parameters
      i=check_lensing_parameters(a);
      if(i)
	{
	  if(!JPL_FLAG)
	    printf("PARAM CHECK: %d %e\n",i,a[abs(i)]);
	  pcheck_count=pcheck_count+1;
	  if(pcheck_count > 1e5){
	    printf("TOO MANY PARAM CHECK (over a million in a row)");  // or else output files can get too big
	    exit(0);
	  }
	  continue;
	}
      // only have these be checks in a row
      pcheck_count = 0;

      // Alexie added test on high mass slope here
      /*
      mh_test=ms_to_mhalo(pow(10.0,11.5), a);
      if(mh_test>1.0E+16){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      mh_test=ms_to_mhalo(pow(10.0,12.0), a);
      if(mh_test>1.0E+20){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      // test the REDs as well
      mh_test=ms_to_mhalo(pow(10.0,11.5),&(a[10]));
      if(mh_test>1.0E+16){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      mh_test=ms_to_mhalo(pow(10.0,12.0), &(a[10]));
      if(mh_test>1.0E+20){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      */

      //xi_2h_gm(1.0,1.0e9,2.0e9,a);
      //chi2wp = chi2smf = chi2lensing = drand48();

      // testing zbrent error
      if(!JPL_FLAG)
	{
	  printf("TEST %d ",count+1);
	  for(i=1;i<=n;++i)
	    printf("%.4e ",a[i]);
	  printf("\n");
	}

      // Calculate chi squares
      chi2smf = chi2_stellar_mass_function(a);

      if(chi2smf<0.99e+07)
	{
	  chi2wp = chi2_wp_lensing_wrapper(a);
	  chi2lensing = chi2_lensing(a);
	  chi2fred = chi2_red_central_fraction(a);
	  chi2 = chi2wp + chi2lensing + chi2smf + chi2fred;
	}
      else
	chi2 = chi2smf;
	 
      printf("\n>> CHI SQR : %e %e %e %e\n",chi2smf,chi2wp,chi2lensing,chi2fred);fflush(stdout);

      if(!ThisTask){
	printf("TRY %d ",++count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("%e %e %e %e\n",chi2smf,chi2wp,chi2lensing,chi2fred);fflush(stdout);
      }

      pcheck[pcnt]=1;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  /* This for loop puts the prev element in the chain is
	   * the current trial point is rejected.
	   */
	  /* For the initialization, don't use this: we need
	   * separate elements for estimating the covariance matrix.
	   */
	  /*
	  for(i=1;i<=n;++i)
	    a[i] = aprev[i];
	  chi2 = chi2prev;
	  */
	  if(USE_IWEIGHT)
	    iweight[niter+1]++;
	  pcheck[pcnt]=0;
	  continue;
	}

      niter++;
      iweight[niter]++;

      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      if(!ThisTask) {
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e %e %e %e %e\n",chi2smf,chi2wp,chi2lensing,chi2fred,chi2);fflush(stdout);
      }

    }

 SKIP_BURN:

  /* hardwire stepfac This is after the first NSTEP runs 1.6/sqrt(n)
   * is about the right value but can change (n=num free parameters)
   * here if the acceptance rate is not 10% < acceptance < 45%
   * optimally, should be about 33%. For a random walk Metropolis,
   * high acceptance rate means that most new samples occur right
   * around the current data point. Their frequent acceptance means
   * that the Markov chain is moving rather slowly and not exploring
   * the parameter space fully. On the other hand, a low acceptance
   * rate means that the proposed samples are often rejected; hence
   * the chain is not moving much. An efficient Metropolis sampler has
   * an acceptance rate that is neither too high nor too low. Roberts
   * and Rosenthal (2001) empirically demonstrated that an acceptance
   * rate between 0.15 and 0.5 is at least 80% efficient, so there is
   * really no need to fine-tune the algorithms to reach acceptance
   * probability that is within small tolerance of the optimal values.
   */
  stepfac=1.6/sqrt(n);
  t0 = second();

  // use an input stepfac if specific
  if(wpl.stepfac>0)
    stepfac *= wpl.stepfac;

  NSTEP = niter;
  pcheck_count = 0;

  /* this is the full loop where the covariance matrix of the
   * chain is used to select the parameters.
   */
  while(niter<imax_chain)
    {

      if(convergence)goto SKIP_MATRIX;
      // if(niter>NSTEP_MAX && niter%NSTEP_MAX!=0)goto SKIP_MATRIX;

      for(j=1;j<=n;++j)
	{
	  avg1[j]=0;
	  for(k=1;k<=n;++k)
	    cov1[j][k]=0;
	}
      total_weight = 0;
      for(i=1;i<=niter;++i)
	{
	  for(j=1;j<=n;++j)
	    {
	      avg1[j]+=chain[i][j]*iweight[i];
	      for(k=1;k<=n;++k)
		cov1[j][k]+=chain[i][j]*chain[i][k]*iweight[i];
	    }
	  total_weight+=iweight[i];
	}

      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp[i][j] = cov1[i][j]/total_weight - avg1[i]*avg1[j]/(total_weight*total_weight);

      jacobi(tmp,n,eval,evect,&nrot);
      gaussj(evect,n,tmp1,1);

    SKIP_MATRIX:
      if(RESTART==4)convergence = 1;

      if(RESTART && count==0)stepfac=0;
      for(i=1;i<=n;++i)
	atemp[i] = gasdev(&IDUM)*sqrt(eval[i])*stepfac;

      for(i=1;i<=n;++i)
	for(a[i]=0,j=1;j<=n;++j)
	  a[i] += atemp[j]*evect[j][i];

      for(i=1;i<=n;++i) 
	a[i] += aprev[i];

      for(i=1;i<=n;++i) 
	wpl.a[i] = a[i];
      wpl.reset_inversion = 1;

      // Check parameters
      i=check_lensing_parameters(a);
      if(i)
	{
	  if(!JPL_FLAG)
	    printf("PARAM CHECK: %d %e\n",i,a[abs(i)]);
	  pcheck_count=pcheck_count+1;
	  if(pcheck_count > 1e5){
	    printf("TOO MANY PARAM CHECK (over a million in a row)");  
	    // or else output files can get too big
	    exit(0);
	  }
	  continue;
	}
      // only have these be checks in a row
      if(pcheck_count)printf("CHECK COUNT: %d\n",pcheck_count);
      pcheck_count = 0;

      /*
      // Alexie added test on high mass slope here
      mh_test=ms_to_mhalo(pow(10.0,11.5), a);
      if(mh_test>pow(10.0,16)){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      mh_test=ms_to_mhalo(pow(10.0,12.0), a);
      if(mh_test>pow(10.0,20)){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      // test the REDs as well
      mh_test=ms_to_mhalo(pow(10.0,11.5),&(a[10]));
      if(mh_test>1.0E+16){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      mh_test=ms_to_mhalo(pow(10.0,12.0), &(a[10]));
      if(mh_test>1.0E+20){
	printf(">> High mass slope is too high, skipping this parameter set ..\n");
	continue;
      }
      */
      //xi_2h_gm(1.0,1.0e9,2.0e9,a);
      //chi2wp = chi2smf = chi2lensing = drand48();

      // Testing this set of parameters
      if(!JPL_FLAG)
	{
	  printf("TEST %d ",count+1);
	  for(i=1;i<=n;++i)
	    printf("%.4e ",a[i]);
	  printf("\n");
	}

      chi2smf = chi2_stellar_mass_function(a);
      if(chi2smf<0.99e+07)
	{
	  chi2wp = chi2_wp_lensing_wrapper(a);
	  chi2lensing = chi2_lensing(a);
	  chi2fred = chi2_red_central_fraction(a);
	  chi2 = chi2wp + chi2lensing + chi2smf + chi2fred;
	}
      else
	chi2 = chi2smf;

      tprev = t0;
      t0 = second();
      printf("TRY %d ",++count);
      for(i=1;i<=n;++i)
	printf("%.4e ",a[i]);
      printf("%e %e %e %e\n",chi2smf,chi2wp,chi2lensing,chi2fred);fflush(stdout);
      
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  /*
	  for(i=1;i<=n;++i)
	    a[i] = aprev[i];
	  chi2 = chi2prev;
	  */
	  if(USE_IWEIGHT)
	    iweight[niter+1]++;
	  continue;
	}

      //      if(NSTEP<NSTEP_MAX)NSTEP++;
      niter++;
      if(!convergence)NSTEP = niter;
      iweight[niter]++;

      if(niter%NSTEP_MAX==0 && !convergence && niter>NSTEP_MAX)
	{
	  convergence = 1;
	  for(i=1;i<=n;++i)
	    {
	      x1=fabs(eval[i]-eval_prev[i])/eval_prev[i];
	      if(x1>0.01)convergence = 0;
	      printf("CONVERGENCE CHECK %d %d %e %e %e\n",niter/NSTEP_MAX,i,x1,eval[i],eval_prev[i]);
	    }
	  for(i=1;i<=n;++i)
	    eval_prev[i] = eval[i];
	  convergence = 0;

	  if(convergence)
	    printf("CONVERGENCE ACCOMPLISHED %d %d \n",niter,count);	    
	}
      if(niter==NSTEP_MAX)
	{
	  for(i=1;i<=n;++i)
	    eval_prev[i] = eval[i];
	}


      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      if(!ThisTask) {
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e %e %e %e %e\n",chi2smf,chi2wp,chi2lensing,chi2fred,chi2);fflush(stdout);
      }
    }
}

/* Read in and store the stellar mass function for this redshift interval.
 * Then calculate the stellar mass function predicted by the model and 
 * return a chi^2
 */
double chi2_stellar_mass_function(double *a)
{
  static int n=0, npr, npb, n_ratio, nr, k, istart_b, istart_r;
  static double *mstar, *nstar, *estar, **covar;
  static double *mstar_b, *nstar_b, *estar_b, **covar_b;
  static double *mstar_r, *nstar_r, *estar_r, **covar_r;
  static double *mstar_ratio, *nstar_ratio, *estar_ratio, **covar_ratio;
  int i, j, i1, ii, ibuf, nlim = 20;
  FILE *fp;
  double dlog10m = 0.2; // Hardwire bin width
  static int niter=0;
  double ng1, ng2, mlo, mhi, ngal, ngalprev, dmsdmh, m1, m2, ns1, ns2, **tmp, **tmp2, mass, nsat;
  char fname[1000];
  double chi2ngal, xmodel[1000], xmodelr[100],xmodelb[100];
  FILE *outfile, *outfile_cov;

  /* The HOD set up is different for the stellar mass function because
   * only really need to know the total number of central and satellites
   * not the full Ncen=f(M) and Nsat=f(sat) (doesn't call these functions)
   */
  if(!n)
    {
      // loop over this two times to get the red/blue data
      // ii==1 -> blue
      // ii==2 -> red
      for(ii=1;ii<=2;++ii)
	{
	  // Alexie changed to local path here (to run at JPL)
	  //sprintf(fname,"/Users/jeremy/SMF/smf_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
	  // JPL
	  //sprintf(fname,"/home/asleauth/HOD/Data/smf_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
	  // Tinker's computer	
	  if(ii==1)
	    sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_SMF/smf_blue_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
	  if(ii==2)
	    sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_SMF/smf_red_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
  
	  if(JPL_FLAG)
	    {
	      if(ii==1)
		sprintf(fname,"./DATA/smf_blue_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
	      if(ii==2)
		sprintf(fname,"./DATA/smf_red_z%.2f_%.2f.dat",wpl.zlo,wpl.zhi);
	    }

	  fp = openfile(fname);
	  n = filesize(fp);
	  
	  mstar = dvector(1,nlim);
	  nstar = dvector(1,nlim);
	  estar = dvector(1,nlim);
	  
	  for(i=1;i<=n;++i)
	    fscanf(fp,"%lf %lf %lf",&mstar[i],&nstar[i],&estar[i]);
	  fclose(fp);
	  fprintf(stderr,"Read [%d] lines from [%s]\n",n,fname);
	  
	  // TEMP TEMP TEMP TEMP	
	  for(i=1;i<=n;++i)
	    nstar[i] = pow(10.0,nstar[i]);
	  
	  // invert the matrix
	  // for now, use same covariance matrix for both (temp)
	  if(COVAR)
	    {	
	      // Now read in the covariance matrix
	      // JPL
	      //sprintf(fname,"/home/asleauth/HOD/Data/smf_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
	      //sprintf(fname,"/Users/alexie/Work/HOD/Data/smf_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
	      // Tinker's computer
	      if(ii==1)
		sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_SMF/smf_blue_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
	      if(ii==2)
		sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_SMF/smf_red_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
  
	      if(JPL_FLAG)
		{
		  if(ii==1)
		    sprintf(fname,"./DATA/smf_blue_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
		  if(ii==2)
		    sprintf(fname,"./DATA/smf_red_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
		}

	      fp = openfile(fname);
	      if(filesize(fp)!=n*n)
		{
		  fprintf(stderr,"FILESIZE mismatch: cij and data not match in SMF: %d %d\n",n*n,filesize(fp));
		  exit(0);
		}
	      covar = dmatrix(1,nlim,1,nlim);
	      for(i=1;i<=n;++i)
		for(j=1;j<=n;++j)
		  fscanf(fp,"%d %d %lf",&i1,&i1,&covar[i][j]);
	      fclose(fp);	
	      fprintf(stderr,"Read [%d] lines from [%s]\n",n*n,fname);

	      //replace the diagonals
	      for(i=1;i<=n;++i)
		estar[i] = sqrt(covar[i][i]);
	      
	      printf("INVERTING COVARIANCE MATRIX\n");
	      tmp=dmatrix(1,n,1,1);
	      tmp2=dmatrix(1,n,1,n);
	      for(i=1;i<=n;++i)
		for(j=1;j<=n;++j)
		  tmp2[i][j]=covar[i][j];
	      gaussj(tmp2,n,tmp,1);
	      for(i=1;i<=n;++i)
		for(j=1;j<=n;++j)
		  covar[i][j]=tmp2[i][j];
	      free_dmatrix(tmp,1,n,1,1);
	      free_dmatrix(tmp2,1,n,1,n);
	      
	      // Write the inverse cov_matrix
	      if(LENSING_OUTPUT_FLAG){
		outfile_cov = fopen("smf_cov_inv.dat","w");
		for(i=1;i<=n;++i)
		  for(j=1;j<=n;++j)
		    fprintf(outfile_cov,"%d %d %e\n",i,j,covar[i][j]);
		fflush(outfile_cov);
		fclose(outfile_cov);
	      }
	    }
	  
	  // now put all the read-in data into the proper data structure
	  if(ii==1) //blue
	    {
	      npb = n;
	      mstar_b = dvector(1,n);
	      nstar_b = dvector(1,n);
	      estar_b = dvector(1,n);
	      covar_b = dmatrix(1,n,1,n);
	      for(i=1;i<=npb;++i)
		{
		  mstar_b[i] = mstar[i];
		  nstar_b[i] = nstar[i];
		  estar_b[i] = estar[i];
		  for(j=1;j<=npb;++j)
		    covar_b[i][j] = covar[i][j];
		}
	    }
	  if(ii==2) //red
	    {
	      npr = n;
	      mstar_r = dvector(1,n);
	      nstar_r = dvector(1,n);
	      estar_r = dvector(1,n);
	      covar_r = dmatrix(1,n,1,n);
	      for(i=1;i<=npr;++i)
		{
		  mstar_r[i] = mstar[i];
		  nstar_r[i] = nstar[i];
		  estar_r[i] = estar[i];
		  for(j=1;j<=npb;++j)
		    covar_r[i][j] = covar[i][j];
		}
	    }
	}
      // now that we have all teh data, calculate the blue/red ratio
      if(wpl.iz==1) 
	{
	  n_ratio = 13;
	  istart_b = 2;
	  istart_r = 1;
	}
      if(wpl.iz==2) 
	{
	  n_ratio = 11;
	  istart_b = 2;
	  istart_r = 1;
	}
      if(wpl.iz==3) 
	{
	  n_ratio = 10;
	  istart_b = 2;
	  istart_r = 1;
	}
      nr = n_ratio;
      mstar_ratio = dvector(1,nr);
      nstar_ratio = dvector(1,nr);
      estar_ratio = dvector(1,nr);
      covar_ratio = dmatrix(1,nr,1,nr);

      for(i=1;i<=n_ratio;++i)
	{
	  j = istart_b+i-1;
	  mstar_ratio[i] = mstar_b[j];
	  nstar_ratio[i] = nstar_b[j];
	  j = istart_r+i-1;
	  nstar_ratio[i] /= nstar_r[j];
	}

      // read in the covariance matrix
      sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_SMF/smf_ratio_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
      if(JPL_FLAG)
	sprintf(fname,"./DATA/smf_ratio_z%.2f_%.2f.covar",wpl.zlo,wpl.zhi);
      
      fp = openfile(fname);
      n = filesize(fp);
      if(n != n_ratio*n_ratio) endrun("ERROR: filesize mismatch for ratio covar!\n");
      n = n_ratio;
      for(i=1;i<=n_ratio;++i)
	for(j=1;j<=n_ratio;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar_ratio[i][j]);
      fclose(fp);	
      fprintf(stderr,"Read [%d] lines from [%s]\n",n_ratio,fname);

      //replace the diagonals
      for(i=1;i<=n;++i)
	estar_ratio[i] = sqrt(covar_ratio[i][i]);

      printf("INVERTING COVARIANCE MATRIX\n");
      tmp=dmatrix(1,n,1,1);
      tmp2=dmatrix(1,n,1,n);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp2[i][j]=covar_ratio[i][j];
      gaussj(tmp2,n,tmp,1);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  covar_ratio[i][j]=tmp2[i][j];
      free_dmatrix(tmp,1,n,1,1);
      free_dmatrix(tmp2,1,n,1,n);
	      

    }

  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.22)<0.1)
    outfile = fopen("smf_outfile.z1","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.48)<0.1)
    outfile = fopen("smf_outfile.z2","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.74)<0.1)
    outfile = fopen("smf_outfile.z3","w");

  chi2ngal = 0;
  ngal = 0;

  for(ii=1;ii<=2;++ii)
    {
      wpl.reset_inversion = 1;
      if(ii==1) //blue
	{
	  BLUE_FLAG = 1;
	  ibuf = 0;
	  n = npb ;
	  for(i=1;i<=npb;++i)
	    {
	      mstar[i] = mstar_b[i];
	      nstar[i] = nstar_b[i];
	      estar[i] = estar_b[i];
	      for(j=1;j<=npb;++j)
		covar[i][j] = covar_b[i][j];
	    }
	}
      if(ii==2) //red
	{
	  BLUE_FLAG = 0;
	  ibuf = HOD_PARAMS;
	  n = npr;
	  for(i=1;i<=npr;++i)
	    {
	      mstar[i] = mstar_r[i];
	      nstar[i] = nstar_r[i];
	      estar[i] = estar_r[i];
	      for(j=1;j<=npr;++j)
		covar[i][j] = covar_r[i][j];
	    }
	}

      i=set_up_hod_for_shmr(1.0E+10+drand48()*1.0E+10,3.0E10,a);
      //if we have an issue, return big chi2 value
      if(i==-1)return 1.0E+7;
	

      for(i=1;i<=n;++i)
	{
	  mlo = pow(10.0,mstar[i]-dlog10m/2); // these are the bin limits (don't use wpl here since Ncen is not called)
	  mhi = pow(10.0,mstar[i]+dlog10m/2);

	  if(set_up_hod_for_shmr(mlo,mhi,a)==-1)return 1.0E+7;
	  ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  nsat = 0;
	  if(LENSING_OUTPUT_FLAG)
	    nsat = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt);

	  // ----
	  
	  ngalprev = ngal;
	  ngal = ngal*pow(HUBBLE,3.0)/dlog10m; // Galaxy density, convert to h72 units
	  ng1 = ng2 = ns1 = ns2 = 0;
      
	  // check for monotonically decreasing mass function
	  if(i>n)
	    {
	      if(ngalprev<ngal)chi2ngal += 1.0E6; // previous density should be higher
	      break;
	    }
	  
	  if(HOD.M_low>HOD.M_max)
	    {
	      if(ERROR_FLAG==1)ERROR_FLAG = 0;
	      if(isnan(ngal))ngal = 0;
	    }
	  
	  if(!COVAR)
	    chi2ngal += (ngal-nstar[i])*(ngal-nstar[i])/(estar[i]*estar[i]);
	  xmodel[i] = ngal;
	  
	  //fprintf(stdout,"SMF%d %e %e %e %e %e %e %e %e %e\n",
	  //	      niter,mstar[i],nstar[i],estar[i],ngal,(ns1-ns2)*pow(HUBBLE,3.0)/dlog10m,ng1,ng2,ns1,ns2);
	  if(LENSING_OUTPUT_FLAG) 
	    {
	      fprintf(outfile,"SMF%d %d %e %e %e %e %e %e %e %e %e\n",
		      niter,ii,mstar[i],nstar[i],estar[i],ngal,nsat*pow(HUBBLE,3.0)/dlog10m,ng1,ng2,ns1,ns2);
	      fflush(outfile);
	    }
	  //  fprintf(stdout,"SMF%d %e %e %e %e %e %e %e %e %e\n",
	  //	  niter,mstar[i],nstar[i],estar[i],ngal,nsat*pow(HUBBLE,3.0)/dlog10m,ng1,ng2,ns1,ns2);
	  if(ERROR_FLAG == 1){
	    printf("\n STOPPED IN chi2_stellar_mass_function %d\n",i);
	    printf("SMF%d %e %e %e %e %e\n",niter,mstar[i],nstar[i],estar[i],ngal,(ns1-ns2)*pow(HUBBLE,3.0)/dlog10m);
	    printf("%e %e %e %e\n",HOD.M_low, HOD.M_min, HOD.M_max, N_sat(HOD.M_max));
	    ERROR_FLAG=0;
	    //exit(0);
	  }
	  
	}
      
      // now calculate the chi^2 with the covariance matrix
      if(COVAR)
	for(i=1;i<=n;++i)
	  for(j=1;j<=n;++j){
	    chi2ngal += (xmodel[i]-nstar[i])*(xmodel[j]-nstar[j])*covar[i][j];
	  }
      if(!JPL_FLAG)
	printf("CHISMFBLUE %e\n",chi2ngal);

      if(ii==1)
	for(i=1;i<=n;++i)
	  xmodelb[i] = xmodel[i];
      if(ii==2)
	for(i=1;i<=n;++i)
	  xmodelr[i] = xmodel[i];
    }
  if(LENSING_OUTPUT_FLAG)
    fclose(outfile);

  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.22)<0.1)
    outfile = fopen("smfr_outfile.z1","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.48)<0.1)
    outfile = fopen("smfr_outfile.z2","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.74)<0.1)
    outfile = fopen("smfr_outfile.z3","w");


  // Now do the ratio
  for(i=1;i<=n_ratio;++i)
    {
      j = istart_b+i-1;
      k = istart_r+i-1;
      xmodel[i] = xmodelb[j]/xmodelr[k];
      if(LENSING_OUTPUT_FLAG) 
	{
	  fprintf(outfile,"%e %e %e %e\n",
		  mstar[i],nstar_ratio[i],(estar_ratio[i]),xmodel[i]);
	  fflush(outfile);
	}
    }
  if(COVAR)
    for(i=1;i<=n_ratio;++i)
      for(j=1;j<=n_ratio;++j)
	chi2ngal += (xmodel[i]-nstar_ratio[i])*(xmodel[j]-nstar_ratio[j])*covar_ratio[i][j];
  if(LENSING_OUTPUT_FLAG)
    fclose(outfile);
  


  // Set this flag to 1 to avoid fitting the SMF
  // This must go and end and not beginning because otherwise galaxy_density is not calcualted properly.
  // ### why was code failing when galaxy density not calculated ??? 
  if(DONT_FIT_SMF)
    chi2ngal=0;

  if(!JPL_FLAG)
    printf("CHISMF %e\n",chi2ngal);
  
  niter++;
  return chi2ngal;
}


/* Read in the clustering data from the 
 * cosmos fields. Analysis is done one redshift bin at a time,
 * so we read in the stellar-mass thresholds for a given redshift.
 */
void wp_lensing_input(void)
{
  int i,j,n,k,i1, istart, ii, nlim = 20;
  FILE *fp;
  char fname[1000], aa[1000];
  double **covar, **tmp, **tmp2;

  wpl.zlo = wpl.zhi = -1;
  if(REDSHIFT>=0.2 || REDSHIFT<0.5)
    {
      wpl.iz = 1;
      wpl.zlo = 0.22;
      wpl.zhi = 0.48;
      wpl.nsamples = 5;
      istart = 1;
    }
  if(REDSHIFT>=0.5 && REDSHIFT <0.75)
    {
      wpl.iz = 2;
      wpl.zlo = 0.48;
      wpl.zhi = 0.74;
      wpl.nsamples = 6;
      istart = 2;
    }
  if(REDSHIFT>=0.75 && REDSHIFT <1.00)
    {
      wpl.iz = 3;
      wpl.zlo = 0.74;
      wpl.zhi = 1.00;
      wpl.nsamples = 6;
      istart = 3;
    }
  if(wpl.zlo<0)
    endrun("ERROR: no proper redshift specified.\n");


  wpl.mstar_threshold = dvector(1,wpl.nsamples);
  wpl.ndata = ivector(1,wpl.nsamples);
  wpl.ndatar = ivector(1,wpl.nsamples);
  wpl.ndatab = ivector(1,wpl.nsamples);
  wpl.ngal = dvector(1,wpl.nsamples);

  wpl.mstar_threshold[1] = 8.8;
  wpl.mstar_threshold[2] = 9.3;
  wpl.mstar_threshold[3] = 9.8;
  wpl.mstar_threshold[4] = 10.3;
  wpl.mstar_threshold[5] = 10.8;
  wpl.mstar_threshold[6] = 11.1;

  // bins for the BLUE SAMPLE
  wpl.mstar_wplo[1][1] = 8.8;
  wpl.mstar_wplo[1][2] = 9.3;
  wpl.mstar_wplo[1][3] = 9.8;
  wpl.mstar_wplo[1][4] = 10.3;
  wpl.mstar_wplo[1][5] = 10.8;
  wpl.mstar_wplo[1][6] = 11.1;

  wpl.mstar_wphi[1][1] = 9.3;
  wpl.mstar_wphi[1][2] = 9.8;
  wpl.mstar_wphi[1][3] = 10.3;
  wpl.mstar_wphi[1][4] = 10.8;
  wpl.mstar_wphi[1][5] = 11.3;
  wpl.mstar_wphi[1][6] = 11.6;

  // bins for the RED SAMPLE
  wpl.mstar_wplo[2][1] = 8.8;
  wpl.mstar_wplo[2][2] = 9.3;
  wpl.mstar_wplo[2][3] = 9.8;
  wpl.mstar_wplo[2][4] = 10.3;
  wpl.mstar_wplo[2][5] = 10.8;
  wpl.mstar_wplo[2][6] = 11.1;

  wpl.mstar_wphi[2][1] = 9.3;
  wpl.mstar_wphi[2][2] = 9.8;
  wpl.mstar_wphi[2][3] = 10.3;
  wpl.mstar_wphi[2][4] = 10.8;
  wpl.mstar_wphi[2][5] = 11.3;
  wpl.mstar_wphi[2][6] = 11.6;

  // loop over BLUE and RED
  for(ii=1;ii<=2;++ii)
    {
      // read in the wp data
      for(i=istart;i<=wpl.nsamples;++i)
	{
	  // Alexie changed here to read in locally
	  //sprintf(fname,"/Users/jeremy/WTHETA/wtheta_jack.mgt%.1f_z%.2f_%.2f",wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	  // JPL
	  //sprintf(fname,"/home/asleauth/HOD/Data/wtheta_jack.mgt%.1f_z%.2f_%.2f",wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	  //sprintf(fname,"/Users/alexie/Work/HOD/Data/wtheta_jack.mgt%.1f_z%.2f_%.2f",wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	  // Tinker's laptop
	  if(ii==1)
	    sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.sfr_dat",
		    wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	  else
	    sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.quench_dat",
		    wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	  if(JPL_FLAG)
	    {
	      if(ii==1)
		sprintf(fname,"./DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.sfr_dat",
			wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	      else
		sprintf(fname,"./DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.quench_dat",
			wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	    }

	  fp = openfile(fname);
	  n = filesize(fp);
	  wpl.ndata[i] = n;
	  wpl.rdata[i] = dvector(1,nlim);
	  wpl.xdata[i] = dvector(1,nlim);
	  wpl.edata[i] = dvector(1,nlim);
	  wpl.covar[i] = dmatrix(1,nlim,1,nlim);
	  
	  for(j=1;j<=n;++j)
	    {
	      fscanf(fp,"%lf %d %lf %lf",&wpl.rdata[i][j],&i1,&wpl.xdata[i][j],&wpl.edata[i][j]);
	      //printf("%e %e %e\n",wpl.rdata[i][j],wpl.xdata[i][j],wpl.edata[i][j]);
	      fgets(aa,1000,fp);
	    }
	  fprintf(stderr,"Read %d lines from [%s]\n",n,fname);
	  fclose(fp);

	  if(COVAR)
	    {
	      // Now read in the covariance matrix
	      //sprintf(fname,"/Users/jeremy/WTHETA/wtheta_covar.mgt%.1f_z%.2f_%.2f",
	      //	  wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	      // JPL
	      //sprintf(fname,"/home/asleauth/HOD/Data/wtheta_covar.mgt%.1f_z%.2f_%.2f",
	      //	  wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	      //sprintf(fname,"/Users/alexie/Work/HOD/Data/wtheta_covar.mgt%.1f_z%.2f_%.2f",
	      //	  wpl.mstar_threshold[i],wpl.zlo,wpl.zhi);
	      // Tinker's laptop
	      if(ii==1)
		sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.sfr_covar",
			wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	      else
		sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.quench_covar",
			wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
	      if(JPL_FLAG)
		{
		  if(ii==1)
		    sprintf(fname,"./DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.sfr_covar",
			    wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
		  else
		    sprintf(fname,"./DATA/wtheta.m%.1f_%.1f_z%.2f_%.2f.quench_covar",
			    wpl.mstar_wplo[ii][i],wpl.mstar_wphi[ii][i],wpl.zlo,wpl.zhi);
		}

	      
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
	      
	      // for z1bin, 1st mbin, reduce errors by x2
	      if(ii==2 && i==istart && wpl.iz==1)
		{
		  for(j=1;j<=n;++j)
		    for(k=1;k<=n;++k)
		      covar[j][k] /= 4;
		}

	      // replace the diagonals
	      for(j=1;j<=n;++j)
		wpl.edata[i][j] = sqrt(covar[j][j]);
	      

	      printf("INVERTING COVARIANCE MATRIX\n");
	      tmp=dmatrix(1,n,1,1);
	      tmp2=dmatrix(1,n,1,n);
	      for(j=1;j<=n;++j)
		for(k=1;k<=n;++k)
		  tmp2[j][k]=covar[j][k];
	      gaussj(tmp2,n,tmp,1);
	      for(j=1;j<=n;++j)
		for(k=1;k<=n;++k)
		  wpl.covar[i][j][k]=tmp2[j][k];
	      free_dmatrix(tmp,1,n,1,1);
	      free_dmatrix(tmp2,1,n,1,n);
	      free_dmatrix(covar,1,n,1,n);
	    }
	  // put this in the right place
	  if(ii==1)
	    {
	      printf("READING WP: %d %d\n",i,n);
	      wpl.ndatab[i] = n;
	      wpl.rdatab[i] = dvector(1,nlim);
	      wpl.xdatab[i] = dvector(1,nlim);
	      wpl.edatab[i] = dvector(1,nlim);
	      wpl.covarb[i] = dmatrix(1,nlim,1,nlim);

	      for(j=1;j<=n;++j)
		{
		  wpl.rdatab[i][j] = wpl.rdata[i][j];
		  wpl.xdatab[i][j] = wpl.xdata[i][j];
		  wpl.edatab[i][j] = wpl.edata[i][j];
		  for(k=1;k<=n;++k)
		    wpl.covarb[i][j][k] = wpl.covar[i][j][k];
		}
	    }
	  if(ii==2)
	    {
	      wpl.ndatar[i] = n;
	      wpl.rdatar[i] = dvector(1,nlim);
	      wpl.xdatar[i] = dvector(1,nlim);
	      wpl.edatar[i] = dvector(1,nlim);
	      wpl.covarr[i] = dmatrix(1,nlim,1,nlim);

	      for(j=1;j<=n;++j)
		{
		  wpl.rdatar[i][j] = wpl.rdata[i][j];
		  wpl.xdatar[i][j] = wpl.xdata[i][j];
		  wpl.edatar[i][j] = wpl.edata[i][j];
		  for(k=1;k<=n;++k)
		    wpl.covarr[i][j][k] = wpl.covar[i][j][k];
		}
	    }
	}
    }
}

/* This is for ALEXIE:
 * Returns the two-halo part of the galaxy-matter
 * cross correlation function at separation r for stellar mass bin mstar_lo to mstar_hi
 * the parameters given in a[].
 * NB-- must call set_up_hod_for_shmr before calling this function
 * NB-- alexie needs to set EXCLUSION=5 and RESET_FLAG_2H in her code.
 * r= comoving h^-1
 */
double xi_2h_gm(double r, double mlo, double mhi, double *a)
{
  return two_halo_real_space(r);
}

double func_satellite_density_xigm(double m)
{
  m = exp(m);
  return m*dndM_interp(m)*Nsat_xigm(m);
}



// Function for Ncen with Varying scatter ...
// Gaussian with varying scatter ...
// ms is already in log units
// ms is stellar mass is log units
double funct_int_vary_scat(double ms)
{
  double res,scatter,sm_scatter,a_fit,b_fit,c_fit,sum;

  // See study_redshift_errors.pro
  if(REDSHIFT>=0.2 || REDSHIFT<0.5)
    {
      //IDL> print,c_low_z1
      //    3.11660    -0.600545    0.0291694
      //IDL> print,c_hi_z1
      //   3.49321    -0.672634    0.0326140
      a_fit=3.11660;
      b_fit=-0.600545;
      c_fit=0.0291694;
    }
  if(REDSHIFT>=0.5 && REDSHIFT <0.75)
    {
      //IDL> print,c_low_z2
      //  2.64875    -0.492551    0.0230964
      //IDL> print,c_hi_z2
      //  4.91459    -0.931206    0.0443085
      a_fit=2.64875;
      b_fit=-0.492551;
      c_fit=0.0230964;
    }
  if(REDSHIFT>=0.75 && REDSHIFT <1.00)
    {
      //IDL> print,c_low_z3
      //  4.15489    -0.784558    0.0372499
      //IDL> print,c_hi_z3
      //  6.74309     -1.27774    0.0607118
      a_fit=4.15489;
      b_fit=-0.784558;
      c_fit=0.0372499;
    }

  //  aa=-0.015000;
  //  bb=0.21000;  
  //sm_scatter = (aa*ms)+bb;                  //stellar mass dependant scatter (lin in log space)

  sm_scatter=a_fit+(ms*b_fit)+(pow(ms,2)*c_fit);

  sum = pow(sm_scatter,2)+pow(wpl.a[6],2);
  scatter = pow(sum,0.5); //sum in quadrature

  // Gaussian in log space
  res = (1.0/(scatter*pow(2*3.14159,0.5))) * exp( -pow(ms - log10(ms_to_mhalo_inversion(mh_scatter_vary)),2)/(2*pow(scatter,2)));
  
  return(res);
}


/****************************************************
 * HOD functions for N_cen using the CLF formulation
 * and the inversion of the M_h(M_gal) formula
 * This is in THRESHOLDS
 ****************************************************/

/* TO DO: make sure wpl.mstar set EVERY TIME it is updated.
   TO DO: Is there a better way to do the bins?
 */

double Ncen_CLF(double m)
{
  static double mstar_prev=-1;
  static double *mh, *ncen, *zz;
  static int n=0;

  int i, imax, ibuf;
  double mlo, mhi, dlogm, a, ms, max;

  ibuf = (1-BLUE_FLAG)*HOD_PARAMS;

  if(!n)
    {
      n=200;
      mh = dvector(1,n);
      ncen = dvector(1,n);
      zz = dvector(1,n);
    }

  if(mstar_prev != wpl.mstar) // Retabulates on wpl.mstar
    {
      mstar_prev = wpl.mstar;
      wpl.mlow_flag=0;

      mlo = log(1.0E8);
      mhi = log(HOD.M_max);
      dlogm = (mhi-mlo)/(n-1);

      ms = wpl.mstar; // assuming that wpl.mstar is in log10 already!!
      imax =  0 ;
      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-1)*dlogm + mlo);
	  if(wpl.ncf>5){
	    ncen[i] = 0.5*(1 - erf((ms - log10(ms_to_mhalo_inversion(mh[i])))/(ROOT2*wpl.a[6+ibuf])));
	  }
	  else{
	    ncen[i] = 0.5*(1 - erf((ms - log10(ms_to_mhalo_inversion(mh[i])))/(ROOT2*aa6)));
	  }

	  //printf("%e %e %e %e %e\n",mh[i],ncen[i],ms,ms_to_mhalo_inversion(mh[i]),wpl.a[6]);
	  max = red_central_fraction(mh[i],wpl.a);
	  if(BLUE_FLAG) max = 1-max;
	  mh[i] = log(mh[i]);
	  if(ncen[i]>1.0E-3 && !imax)imax = i;
	  ncen[i] = log(ncen[i]*max + 1.0E-20);
	}
      spline(mh,ncen,n,1.0E+30,1.0E+30,zz);
      if(!imax)wpl.mlow_flag=1;
    }

  m = log(m);
  //printf("NCEN %e %e\n",m,mh[1]);
  if(m<mh[1])return 1.0E-20;
  splint(mh,ncen,zz,n,(m),&a);
  if(a>0)return 1;
  return exp(a);
  }



/* loop over the different bins in stellar mass
 * and calculate the correlation function and chi^2
 * for each bin.
 */
double chi2_wp_lensing_wrapper(double *a)
{
  int i,j,k,n,istart,ii,ibuf;
  double dmsdmh,x0,mstar,m1,m2,dlogm,model[100];
  double chi2w, chi2tot, chi2ngal;
  static int niter=0;

  double dlogr, r, rlo, rhi, mlo, mhi;
  FILE *outfile, *fp2;

  if(ERROR_FLAG == 1){
    printf("\n UNCAUGHT ERROR: STOPPED IN chi2_wp_lensing_wrapper\n"); 
  }

  chi2tot = 0;
  istart = 3;
  if(REDSHIFT<0.74)istart = 2;
  if(REDSHIFT<0.48)istart = 1;
  

  if(DONT_FIT_CLUSTERING)
    return 0.0;

  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.22)<0.1)
    outfile = fopen("wtheta_outfile.z1","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.48)<0.1)
    outfile = fopen("wtheta_outfile.z2","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.74)<0.1)
    outfile = fopen("wtheta_outfile.z3","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.22)<0.1)
    fp2 = fopen("hod_outfile.z1","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.48)<0.1)
    fp2 = fopen("hod_outfile.z2","w");
  if(LENSING_OUTPUT_FLAG && fabs(wpl.zlo-0.74)<0.1)
    fp2 = fopen("hod_outfile.z3","w");

  // do a loop for both blue and red
  for(ii=1;ii<=2;++ii)
    {
      if(ii==1)
	{
	  BLUE_FLAG=1;
	  ibuf = 0;
	}
      else
	{
	  BLUE_FLAG=0;
	  ibuf = HOD_PARAMS;
	}
      
      
      for(i=istart;i<=wpl.nsamples;++i)
	{
	  if(i==istart) { wpl.reset_inversion = 1; ms_to_mhalo_inversion(1.0e10); } 
	  mstar = pow(10.0,wpl.mstar_threshold[i]); //### what is this for ?


	  // put the data in place:
	  if(ii==1)
	    {
	      n = wpl.ndata[i] = wpl.ndatab[i];
	      for(j=1;j<=n;++j)
		{
		  wpl.rdata[i][j] = wpl.rdatab[i][j];
		  wpl.xdata[i][j] = wpl.xdatab[i][j];
		  wpl.edata[i][j] = wpl.edatab[i][j];
		  for(k=1;k<=n;++k)
		    wpl.covar[i][j][k] = wpl.covarb[i][j][k];
		  //printf(">> %e %e %e \n",wpl.rdata[i][j], wpl.xdata[i][j], wpl.edata[i][j]);
		}
	    }
	  if(ii==2)
	    {
	      n = wpl.ndata[i] = wpl.ndatar[i];
	      for(j=1;j<=n;++j)
		{
		  wpl.rdata[i][j] = wpl.rdatar[i][j];
		  wpl.xdata[i][j] = wpl.xdatar[i][j];
		  wpl.edata[i][j] = wpl.edatar[i][j];
		  for(k=1;k<=n;++k)
		    wpl.covar[i][j][k] = wpl.covarr[i][j][k];
		}
	    }	


	  mlo = pow(10.0,wpl.mstar_wplo[ii][i]);
	  mhi = pow(10.0,wpl.mstar_wphi[ii][i]);
	  set_up_hod_for_shmr(mlo,mhi,a);

	  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

	  if(ERROR_FLAG == 1){
	    printf("\n STOPPED IN chi2_wp_lensing_wrapper HERE 2\n"); 
	    printf("%e %e %e\n",m1,m2,mstar);
	    printf("FSAT%d %d %f %e %e %e\n",i,niter,qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY,HOD.M_min,HOD.M_low,HOD.M1);
	    exit(0);
	  }
	  
	  RESET_FLAG_1H = RESET_FLAG_2H = 1;
	  
	  chi2w = 0;
	  
	  for(j=1;j<=wpl.ndata[i];++j)
	    {
	      x0 = wtheta(wpl.rdata[i][j])*integral_constraint_correction(wpl.rdata[i][j]);
	      model[j] = x0;


	      // mutiple by bias of 10^13.75 Msol halo
	      mlo = pow(10.0,13.75)*HUBBLE;
	      x0 *= bias_interp(mlo,-1)*bias_interp(mlo,-1);
	      fmuh(bias_interp(mlo,-1));

	      if(LENSING_OUTPUT_FLAG)
		fprintf(outfile,"WTH %d %d %d %e %e %e %e\n",ii,i,niter,wpl.rdata[i][j],wpl.xdata[i][j],x0,wpl.edata[i][j]);
	      if(LENSING_OUTPUT_FLAG)
		fflush(outfile);
	      // fprintf(stdout,"WTH%02d%02d %d %e %e %e %e\n",ii,i,niter,wpl.rdata[i][j],wpl.xdata[i][j],x0,wpl.edata[i][j]);
	      
	      // Alexie commented here to limit ourput	  
	      printf("WTH%d %d %e %e %e %e\n",i,niter,wpl.rdata[i][j],wpl.xdata[i][j],x0,wpl.edata[i][j]);
	      if(!COVAR)
		chi2w += (wpl.xdata[i][j]-x0)*(wpl.xdata[i][j]-x0)/(wpl.edata[i][j]*wpl.edata[i][j]);

	    }
	  // temp for bias
	  exit(0);

	  if(LENSING_OUTPUT_FLAG)
	    {
	      x0 = 1.0e10/1.5;
	      for(j=0;j<1000;++j)
		{
		  x0 = x0*1.2;
		  fprintf(fp2,"HOD %d %d %e %e %e\n",ii,i,x0,N_cen(x0),N_sat(x0));
		  if(x0>1.0E15)break;
		}
	    }

	  
	  if(COVAR)
	    {
	      for(j=1;j<=wpl.ndata[i];++j)
		for(k=1;k<=wpl.ndata[i];++k)
		  chi2w += (model[j]-wpl.xdata[i][j])*(model[k]-wpl.xdata[i][k])
		    *wpl.covar[i][j][k];
	    }
	  
	  //Alexie commented here to limit output
	  if(!JPL_FLAG)
	    {
	      printf("CHIWP%d %d %e\n",i,niter,chi2w);
	      fflush(stdout);
	    }
	  // if this is the first red sample, do not include!
	  //if(!BLUE_FLAG && i==istart)continue; //(nah... keep it!)
	  // if this is the last blue sample for Z2, do not include!
	  if(BLUE_FLAG && i==6 && istart==2){  continue; }
	  chi2tot += chi2w;
	}
    }

  if(LENSING_OUTPUT_FLAG)
    fclose(outfile);
  if(LENSING_OUTPUT_FLAG)
    fclose(fp2);


  printf("CHIWP%d %d %e\n",i,niter,chi2tot);

  niter++;
  return chi2tot;
}


/* check and make sure the parameters selected are within the 
 * acceptable parameter space.
 */
int check_lensing_parameters(double *a)
{
  int check_flag, ibuf;
  float m1;

  wpl.reset_fred = 1;

  if(a[1]<=12.0)return 1;      // Mhalo_norm
  if(a[1]>=13.4)return -1;
    
  if(a[2]<=10.3)return 2;      // Mstellar norm
  if(a[2]>=12.0)return -2;
 
  if(a[3]<=0.2)return 3;       // Beta 
  if(a[3]>=0.8)return -3;

  if(a[4]<=0.1)return 4;       // delta (If <0, affects the threshold HOD, strange things happen .. don't go below 0.3)
  if(a[4]>=1.4)return -4;
  
  if(a[5]<=-0.5)return 5;      // gamma (can be <0)
  if(a[5]>=5.0)return -5;

  if(wpl.ncf>5) {
    if(a[6]<=0.01)return 6;  // sigma_log_sm
    if(a[6]>=0.5)return -6;
  }
  
  if(wpl.ncf>9) {
    
    // This now varies around 2
    if(a[7]<0.0)return 7;     // Norm Mcut (poorly constrained)
    if(a[7]>=40)return -7;

    // This now varies around 17
    if(a[8]<0.0)return 8;     // Norm M1
    if(a[8]>=140)return -8;

    // a[9] (slope) can be negative !!
    if(a[9]<=-2.0)return 9;      // Slope Mcut
    if(a[9]>=6)return -9;

    if(a[10]<=-2.0)return 10;    // Slope M1
    if(a[10]>=4)return -10;
  }

  // 11 parameter model
  if(ALPHA_FREE){
    if(a[11]>2)return 11;
    if(a[11]<0.3)return -11;
  }

  /**-------------------------------------------
   * Now do the parameters for the RED gals
   **------------------------------------------*/
  
  ibuf = HOD_PARAMS;
    
  if(a[ibuf+1]<=11.5)return 1+ibuf;      // Mhalo_norm
  if(a[ibuf+1]>=13.4)return -1-ibuf;
    
  if(a[ibuf+2]<=10.3)return 2+ibuf;      // Mstellar norm
  if(a[ibuf+2]>=12.0)return -2-ibuf;
 
  if(a[ibuf+3]<=0.01)return 3+ibuf;       // Beta 
  if(a[ibuf+3]>=0.8)return -3-ibuf;

  if(a[ibuf+4]<=0.1)return 4+ibuf;       // delta (If <0, affects the threshold HOD, strange things happen .. don't go below 0.3)
  if(a[ibuf+4]>=1.4)return -4-ibuf;
  
  if(a[ibuf+5]<=-0.5)return 5+ibuf;      // gamma (can be <0)
  if(a[ibuf+5]>=5.0)return -5-ibuf;

  if(wpl.ncf>5) {
    if(a[ibuf+6]<=0.01)return 6+ibuf;  // sigma_log_sm
    if(a[ibuf+6]>=0.5)return -6-ibuf;
  }

  if(wpl.ncf>9) {

    // This now varies around 2
    if(a[ibuf+7]<0.0)return 7+ibuf;     // Norm Mcut (poorly constrained)
    if(a[ibuf+7]>=40)return -7-ibuf;

    // This now varies around 17
    if(a[ibuf+8]<0.0)return 8+ibuf;     // Norm M1
    if(a[ibuf+8]>=140)return -8-ibuf;

    // a[ibuf+9] (slope) can be negative !!
    if(a[ibuf+9]<=-2.0)return 9+ibuf;      // Slope Mcut
    if(a[ibuf+9]>=6)return -9-ibuf;

    if(a[ibuf+10]<=-2.0)return 10+ibuf;    // Slope M1
    if(a[ibuf+10]>=4)return -10-ibuf;
  }

  // 11 parameter model
  if(ALPHA_FREE){
    if(a[ibuf+11]>2)return 11+ibuf;
    if(a[ibuf+11]<0.3)return -11-ibuf;
  }

  if(a[23]>a[24])return -23;
  if(a[24]>a[25])return -24;
  if(a[25]>a[26])return -25;
  if(a[26]>a[27])return -26;
  if(a[23]<-10)return 23;
  if(a[24]<-10)return 24;

  return 0;
}

/* set the initial values of the paramters in the 
 * a vector.
 * Set the initial dispersions for burn-in parameter selection.
 */
double mcmc_lensing_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0, i1, i2, n;
  double x1,x2,x3,x4,omega_m, t0, t1;
  long IDUM = -556;

  FILE *infile, *fp;
  char aa[1000];

  double rlo, rhi, dlogr, r,mh_test, mass;
  double *eval, **evect;

  eval = dvector(1,wpl.ncf);
  evect = dmatrix(1,wpl.ncf, 1,wpl.ncf);

  wpl.zlo = wpl.zhi = -1;
  if(REDSHIFT>=0.2 || REDSHIFT<0.5)
    {
      wpl.zlo = 0.22;
      wpl.zhi = 0.48;
    }
  if(REDSHIFT>=0.5 && REDSHIFT <0.75)
    {
      wpl.zlo = 0.48;
      wpl.zhi = 0.74;
    }
  if(REDSHIFT>=0.75 && REDSHIFT <1.00)
    {
      wpl.zlo = 0.74;
      wpl.zhi = 1.00;
    }
  if(wpl.zlo<0)
    endrun("ERROR: no proper redshift specified.\n");



  /*;------------------------------------------------------------------------
    ;  mhalo_norm     = a[1]      
    ;  mstellar_norm  = a[2]   
    ;  beta           = a[3]       beta
    ;  delta          = a[4]       delta
    ;  gamma          = a[5]       gamma
    ;  sigma_logsm    = a[6]       sigma
    ;  Mcut norm      = a[7]       Mcut norm
    ;  M1 norm        = a[8]       M1 norm
    ;  Mcut slope     = a[9]       Mcut slope 
    ;  M1 slope       = a[10]      M1 slope
    ;  sat slope      = a[11]
    ;------------------------------------------------------------------------*/

  // For the 6p model
  aa6  = 0.2;
  aa7  = 2;
  aa8  = 17;
  aa9  = 1.3;
  aa10 = 0.8;
  // Default: set the slope to 1
  HOD.alpha = 1;

  if(REDSHIFT>=0.2 && REDSHIFT<0.5)
    {
      start_dev[1]=0.04;                   // Mhalo norm
      start_dev[2]=0.03;                   // Mstellar norm
      start_dev[3]=0.043;                   // beta
      start_dev[4]=0.10;                   // delta
      start_dev[5]=0.2;                    // gamma
      if(wpl.ncf>5) start_dev[6]=0.05;     // sigma 
      if(wpl.ncf>9){                       // 10p model
	start_dev[7] = 0.9;                // Mcut norm
	start_dev[8] = 3;                  // M1 norm
	start_dev[9] = 0.9;                // Mcut slope
	start_dev[10] = 0.1; }            // M1 slope
      if(ALPHA_FREE) start_dev[11] = 0.05; // 11p model: vary the satellite slope
      
      a[1] = 12.5;                      // Mhalo norm
      a[2] = 10.9;                      // Mstelar norm
      a[3] = 0.46;                      // beta
      a[4] = 0.47;                      // delta
      a[5] = 1.6;                      // gamma
      if(wpl.ncf>5) a[6] = 0.18;        // Sigma
      if(wpl.ncf>9){                    // 10 parameter model
	a[7] = 0.1;                     // Mcut norm
	a[8] = 15.0;                      // M1 norm
	a[9] = 2.0;                    // Mcut slope 
	a[10] = 0.6;}                   // M1 slope
      if(ALPHA_FREE) a[11] = 1.0 ;
    }

  if(REDSHIFT>=0.5 && REDSHIFT <0.75)
    {
      start_dev[1]=0.04;                   // Mhalo norm
      start_dev[2]=0.03;                   // Mstellar norm
      start_dev[3]=0.02;                   // beta
      start_dev[4]=0.12;                   // delta
      start_dev[5]=0.25;                    // gamma
      if(wpl.ncf>5) start_dev[6]=0.02;     // sigma 
      if(wpl.ncf>9){                       // 10p model
	start_dev[7] = 1.3;                // Mcut norm
	start_dev[8] = 2.0;                  // M1 norm
	start_dev[9] = 0.2;                // Mcut slope
	start_dev[10] = 0.08; }            // M1 slope
      if(wpl.ncf>10) start_dev[11] = 0.05; // 11p model: vary the satellite slope
      
      a[1] = 12.704;                     // Mhalo norm
      a[2] = 11.017;                     // Mstelar norm
      a[3] = 0.467;                      // beta
      a[4] = 0.6;                        // delta
      a[5] = 1.9;                        // gamma
      if(wpl.ncf>5) a[6] = 0.247;        // Sigma
      if(wpl.ncf>9){                     // 10 parameter model
	a[7] = 2.0;                      // Mcut norm
	a[8] = 9.0;                        // M1 norm
	a[9] = 0.9;                      // Mcut slope 
	a[10] = 0.7;}                    // M1 slope
      if(wpl.ncf>10) a[11] = 1.0 ;
    }

  if(REDSHIFT>=0.75 && REDSHIFT <1.00)
    {
      start_dev[1]=0.03;                   // Mhalo norm
      start_dev[2]=0.03;                   // Mstellar norm
      start_dev[3]=0.03;                   // beta
      start_dev[4]=0.09;                   // delta
      start_dev[5]=0.3;                    // gamma
      if(wpl.ncf>5) start_dev[6]=0.03;     // sigma 
      if(wpl.ncf>9){                       // 10p model
	start_dev[7] = 0.5;                // Mcut norm
	start_dev[8] = 2.0;                // M1 norm
	start_dev[9] = 0.3;                // Mcut slope
	start_dev[10] = 0.2; }             // M1 slope
      if(wpl.ncf>10) start_dev[11] = 0.05; // 11p model: vary the satellite slope
      
      a[1] = 12.70;                       // Mhalo norm
      a[2] = 11.140;                      // Mstelar norm
      a[3] = 0.4602;                      // beta
      a[4] = 0.39;                        // delta
      a[5] = 2.5;                         // gamma
      if(wpl.ncf>5) a[6] = 0.23;          // Sigma
      if(wpl.ncf>9){                      // 10 parameter model
	a[7] = 2.0 ;                      // Mcut norm
	a[8] = 9 ;                        // M1 norm
	a[9] = 0.5;                       // Mcut slope 
	a[10] = 0.87 ;}                   // M1 slope
      if(wpl.ncf>10) a[11] = 1.0 ;
    }

  //-------------- END OF INITIALIZE VALUES -----
  
  if(LENSING_OUTPUT_FLAG==-1)
    {
      if(REDSHIFT<=0.48)
	infile = fopen("/Users/alexie/Work/HOD/Minchi2/z1_minchi2.txt","r");
      if(REDSHIFT>0.48 && REDSHIFT <=0.74)
	infile = fopen("/Users/alexie/Work/HOD/Minchi2/z2_minchi2.txt","r");
      if(REDSHIFT>0.74 && REDSHIFT <=1.00)
	infile = fopen("/Users/alexie/Work/HOD/Minchi2/z3_minchi2.txt","r");

      fscanf(infile,"%s %d %d",aa,&i1,&i2);
      for(i=1;i<=wpl.ncf;++i)
	fscanf(infile,"%lf",&a[i]);
    }


  // Test model (must put here)
  // Paper FITS FOR z1 :
  a[1] = 12.620;                     // Mhalo norm
  a[2] = 10.916;                     // Mstelar norm
  a[3] = 0.457;                      // beta
  a[4] = 0.766;                      // delta
  a[5] = 1.53;                       // gamma
  a[6] = 0.206;                      // Sigma
  a[7] = 1.47;                       // Mcut norm
  a[8] = 20.62 ;                     // M1 norm
  a[9] = -0.13;                      // Mcut slope 
  a[10] = 0.859 ;                    // M1 slope
  if(ALPHA_FREE)a[11] = 0.8;         // alpha_sat


  // SHMR for the reds
  a[11+ALPHA_FREE] = 12.220;                     // Mhalo norm
  a[12+ALPHA_FREE] = 10.85;                     // Mstelar norm
  a[13+ALPHA_FREE] = 0.357;                      // beta
  a[14+ALPHA_FREE] = 0.766;                      // delta
  a[15+ALPHA_FREE] = 1.53;                       // gamma
  a[16+ALPHA_FREE] = 0.206;                      // sigma
  a[17+ALPHA_FREE] = 1.47;                       // Mcut norm
  a[18+ALPHA_FREE] = 12.62 ;                     // M1 norm
  a[19+ALPHA_FREE] = -0.13;                      // Mcut slope 
  a[20+ALPHA_FREE] = 0.7 ;                    // M1 slope
  if(ALPHA_FREE)a[22]= 1.2;          // alpha_sat

  start_dev[11+ALPHA_FREE] = start_dev[1];
  start_dev[12+ALPHA_FREE] = start_dev[2];
  start_dev[13+ALPHA_FREE] = start_dev[3];
  start_dev[14+ALPHA_FREE] = start_dev[4];
  start_dev[15+ALPHA_FREE] = start_dev[5];
  start_dev[16+ALPHA_FREE] = start_dev[6];
  start_dev[17+ALPHA_FREE] = start_dev[7];
  start_dev[18+ALPHA_FREE] = start_dev[8];
  start_dev[19+ALPHA_FREE] = start_dev[9];
  start_dev[20+ALPHA_FREE] = start_dev[10];
  if(ALPHA_FREE)start_dev[22] = start_dev[11];

  // red central fraction
  a[21+ALPHA_FREE*2] = -2;
  a[22+ALPHA_FREE*2] = -1.5;
  a[23+ALPHA_FREE*2] = 0.4;
  a[24+ALPHA_FREE*2] = 0.5;
  a[25+ALPHA_FREE*2] = 0.6;

  for(i=21+ALPHA_FREE*2;i<=25+ALPHA_FREE*2;++i)
    start_dev[i] = a[i]*0.2;

  // other parameters (not needed right now).

  // Test Behroozi model
  // Note: STILL NEED TO CONVERT MVIR HERE ...
  /*--- Behroozi, z=0.1 -
    ms       10.6700
    mh       12.3245
    beta      0.423636
    delta      0.554545
    gamma       1.33182
    ---------------------*/

  if(LENSING_OUTPUT_FLAG || MAKE_MOCKS)
    {
      fp = openfile(ARGV[2]);
      fscanf(fp,"%s %d %d",aa,&i,&j);
      for(i=1;i<=wpl.ncf;++i)
	fscanf(fp,"%lf",&a[i]);
      fclose(fp);
    }

  if(ARGC>3)
    {
      fp = openfile(ARGV[3]);
      fscanf(fp,"%s %d %d",aa,&i,&j);
      for(i=1;i<=wpl.ncf;++i)
	fscanf(fp,"%lf",&a[i]);
      fclose(fp);
    }

  // if a[23] is <-10 (new contraint) then up it to -9+/-0.5
  if(a[23]<-10)
    a[23] = -8.5-drand48();


  for(i=1;i<=wpl.ncf;++i)
    wpl.a[i] = a[i];

  // this is the test for L120b
  /*
  sprintf(Files.HaloFile,"../L120b_MOCK/halo_fof.2_083.dat");
  BLUE_FLAG = 0;
  wpl.reset_inversion = 1;
  set_up_hod_for_shmr(1.0E+8*drand48(),2.0E+8, a);
  sprintf(Task.root_filename,"../L120b_MOCK/z1_red_fof.2",i);	    
  populate_simulation_old();
  exit(0);
  */  

  /* flag for shuffling things around a bit to get 
   * the starting parameters away for minimum
   */
  if(FISHER==3)
    fisher_matrix(a,eval,evect,1);


  if(MAKE_MOCKS)
    {
      for(i=1;i<=wpl.ncf;++i)
	printf("a[%02d] = %e\n",i,wpl.a[i]);

      
      // output the HOD for reds 9 to 9.5
      wpl.reset_inversion = 1;
      BLUE_FLAG = 0;
      set_up_hod_for_shmr(pow(10.0,11.0),pow(10.0,11.5),a);
      for(i=110;i<=150;++i)
	{
	  mass = pow(10.0,i/10.0);
	  printf("RED %e %e %e\n",mass,N_cen(mass),N_sat(mass));
	}
      printf("BOO %e\n",HOD.M_min);

      // reset the tabulation of Ncen_xigm
      set_up_hod_for_shmr(1.0E+8,2.0E+8, a);

      // output the HOD for blues 9 to 9.5
      wpl.reset_inversion = 1;
      BLUE_FLAG = 1;
      set_up_hod_for_shmr(pow(10.0,11.0),pow(10.0,11.5),a);
      for(i=110;i<=150;++i)
	{
	  mass = pow(10.0,i/10.0);
	  printf("BLUE %e %e %e\n",mass,N_cen(mass),N_sat(mass));
	}
      exit(0);
      

      /* MOCKS MOCKS MOCKS
       */
      // make the mocks!!!
      // 405 for z1
      // 172 for z2
      // 109 for z3
      for(i=0;i<405;++i)
	{
	  sprintf(Files.HaloFile,
		  "/Users/tinker/cosmo/COSMOS/MOCKS/FOF/survey_z0.24-0.48_x91.00_y91.00_%d.dat",i);
	  //"/Users/jeremy/MOCKS/FOF/survey_z0.24-0.48_x91.00_y91.00_%d.dat",i);
	  //"/Users/jeremy/MOCKS/FOF/survey_z0.48-0.74_x91.00_y91.00_%d.dat",i);
	  //"/Users/jeremy/MOCKS/FOF/survey_z0.74-1.00_x91.00_y91.00_%d.dat",i);
	  BLUE_FLAG = 0;
	  wpl.reset_inversion = 1;
	  set_up_hod_for_shmr(1.0E+8*drand48(),2.0E+8, a);
	  sprintf(Task.root_filename,"/Users/tinker/cosmo/COSMOS/MOCKS/GALMOCKS_QUENCHED/z1_red.%d",i);	    
	  //sprintf(Task.root_filename,"/Users/jeremy/MOCKS/GALMOCKS_QUENCHED/z1_red.%d",i);
	  populate_simulation();


	  BLUE_FLAG = 1;
	  wpl.reset_inversion = 1;
	  set_up_hod_for_shmr(1.0E+8*drand48(),2.0E+8, a);
	  sprintf(Task.root_filename,"/Users/tinker/cosmo/COSMOS/MOCKS/GALMOCKS_QUENCHED/z1_blue.%d",i);
	  //sprintf(Task.root_filename,"/Users/jeremy/MOCKS/GALMOCKS_QUENCHED/z1_blue.%d",i);
	  populate_simulation();
	}
      exit(0);

      for(i=0;i<10;++i)
	{
	  sprintf(Files.HaloFile,
		  "/Users/jeremy/MOCKS/LARGE_FOF/survey_z0.00-1.20_x77.00_y77.00_%d.dat",i);
	  sprintf(Task.root_filename,"/Users/jeremy/MOCKS/GALMOCKS/lightcone.%d",i);
	  populate_simulation();
	}
      exit(0);
    }




  
  printf("\n>> LENSING OUTPUT FLAG =  %i \n",LENSING_OUTPUT_FLAG);
  printf(">> DONT FIT LENSING FLAG =  %i \n",DONT_FIT_LENSING);
  printf(">> DONT FIT CLUSTERING FLAG =  %i \n",DONT_FIT_CLUSTERING);
  printf(">> DONT FIT SMF FLAG =  %i \n",DONT_FIT_SMF);
  printf(">> REDSHIFT =  %f \n",REDSHIFT);
  printf(">> COVAR =  %i \n\n",COVAR);

  printf("THE NUMBER OF PARAMS IS :  ");
  muh(wpl.ncf);
  printf("\n");
  for(i=1;i<=wpl.ncf;++i)
    printf(">> a[%02d] = %e\n",i,wpl.a[i]);
  
  printf("ACCEPT 1 1 ");
  for(i=1;i<=wpl.ncf;++i)
    printf("%e ",wpl.a[i]);
  printf(" 0.0 0.0 0.0 0.0\n");

  // Changed to 11.5 here
  mh_test=ms_to_mhalo(pow(10.0,11.5), a);
  printf("\n >>> Mh_TEST at 10^11.5 %e :\n",mh_test);


  // chain_post_process_clf(a);

  if(LENSING_OUTPUT_FLAG)
    chain_post_process_single(a);

  // SMF
  printf("\n >>> Reading SMF for first time :\n");
  x3 = chi2_stellar_mass_function(a);

  // For the wp
  t0 = second();
  printf(" >>> Reading WP for first time :\n");
  x1 = chi2_wp_lensing_wrapper(a);
  t1 = second();
  printf("TIMING: %.2f\n",timediff(t0,t1));//exit(0);

  // ALEXIE: CHI FOR LENSING IS CALLED HERE
  wpl.reset_inversion = 1;
  printf(" >>> Reading LENSING for first time :\n");

  t0 = second();
  x2 = chi2_lensing(a);
  t1 = second();
  printf("TIMING: %.2f\n",timediff(t0,t1));

  wpl.reset_inversion = 1;
  t0 = second();
  x2 = chi2_lensing(a);
  t1 = second();
  printf("TIMING: %.2f\n",timediff(t0,t1));

  x4 = chi2_red_central_fraction(a);


  if(LENSING_OUTPUT_FLAG)
    exit(0);

  printf("\n >>> INITIAL CHI2 (smf, wp, lensing): %e %e %e %e\n",x3,x1,x2,x4);
  printf(" >>> INITIAL CHI2 TOTAL: %e \n\n",x3+x1+x2+x4);

  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wpl.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e %e %e %e\n",x3,x1,x2,x4);fflush(stdout);
    printf("INITIAL CHI2: %e %e %e %e\n",x3,x1,x2,x4);
    fflush(stdout);
  }
  return(x1+x2+x3+x4);
}


double integral_constraint_correction(double r)
{
  //  if($iz==1) { set x = exp((r/3.4)**2.5) } ## for z1
  //if($iz==2) { set x = exp((r/3.2)**4.1) } ## for z2 (alternate)
  //if($iz==3) { set x = r*0 + 1 } ## for z3

  r = log10(r);

  if(REDSHIFT>0.1 && REDSHIFT<0.4)
    return 1.0/exp(pow(r/3.4,2.5));
  if(REDSHIFT>0.5 && REDSHIFT<0.7)
    return 1.0/exp(pow(r/3.2,4.1));
  if(REDSHIFT>0.75)
    return 1;

  printf("ERROR: REDSHIFT=%f not in range for IC correction\n",REDSHIFT);
  exit(0);

  return 0;
}

/**********************************************************************
 * DATA: quenched central fraction
 **********************************************************************/

double chi2_red_central_fraction(double *a)
{
  static int flag=1, nhalo, niter=0;
  static double fqdata, eqdata, *mass, *redshift;

  FILE *fp;
  int i,j,n;
  double x1, chi2, ftot;
  
  if(flag)
    {
      flag = 0;

      // first, set the overall quenched fraction and the error
      if(fabs(wpl.zlo-0.22)<0.1) { fqdata = 0.641; eqdata = 0.056; }
      if(fabs(wpl.zlo-0.44)<0.1) { fqdata = 0.709; eqdata = 0.060; }
      if(fabs(wpl.zlo-0.74)<0.1) { fqdata = 0.643; eqdata = 0.069; }

      fp = openfile("central_galaxies.dat");
      nhalo = filesize(fp);
      mass = dvector(1,nhalo);
      redshift = dvector(1,nhalo);

      for(i=1;i<=nhalo;++i)
	fscanf(fp,"%lf %lf %lf %lf %lf %d",&mass[i],&redshift[i],&x1,&x1,&x1,&j);
      fclose(fp);

      for(i=1;i<=nhalo;++i)
	mass[i] = pow(10.0,mass[i])*HUBBLE; //conver to h-inverse units
    }

  ftot = 0;
  n = 0;
  for(i=1;i<=nhalo;++i)
    {
      if(redshift[i]>wpl.zhi || redshift[i]<wpl.zlo)continue;
      n++;
      ftot += red_central_fraction(mass[i],a);
      //printf("%e %e\n",mass[i],red_central_fraction(mass[i],a));
    }
  ftot /= n;
  chi2 = (ftot-fqdata)*(ftot-fqdata)/(eqdata*eqdata);
  printf("CHIFCEN %d %e %e %e %e\n",niter++,chi2,ftot,fqdata,eqdata);
  return chi2;

}

/*-----------------------------------------------------------------------------
 * 
 * INPUT AN EXTERNAL CHAIN FOR COVARIANCE MATRIX
 *
 *----------------------------------------------------------------------------*/

void input_chain(double *eval, double **evect)
{
  int np,i,j,n,total_weight,nrot,i1,k;
  double **chain, **cov1, *avg1, **tmp, **tmp1;
  FILE *fp;
  char aa[1000];

  fp = openfile(ARGV[4]);
  np = filesize(fp);
  chain = dmatrix(1,np,1,wpl.ncf);

  n = wpl.ncf;
  cov1=dmatrix(1,n,1,n);
  avg1=dvector(1,n);

  tmp=dmatrix(1,n,1,n);
  tmp1=dmatrix(1,n,1,1);

  muh(0);
  for(i=1;i<=np;++i)
    {
      fscanf(fp,"%s %d %d",aa,&i1,&i1);
      for(j=1;j<=wpl.ncf;++j)
	fscanf(fp,"%lf",&chain[i][j]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  muh(1);

  for(j=1;j<=n;++j)
    {
      avg1[j]=0;
      for(k=1;k<=n;++k)
	cov1[j][k]=0;
    }
  total_weight = 0;
  for(i=1;i<=np;++i)
    {
      for(j=1;j<=n;++j)
	{
	  avg1[j]+=chain[i][j];
	  for(k=1;k<=n;++k)
	    cov1[j][k]+=chain[i][j]*chain[i][k];
	}
      total_weight++;
    }
  
  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      tmp[i][j] = cov1[i][j]/total_weight - avg1[i]*avg1[j]/(total_weight*total_weight);
  
  jacobi(tmp,n,eval,evect,&nrot);
  gaussj(evect,n,tmp1,1);

  for(i=1;i<=n;++i)
    printf("EVAL %d %e\n",i,eval[i]);
  
}



/*-----------------------------------------------------------------------------
 * 
 * FISHER MATRIX CALCULATION OF THE COVARIANCE MATRIX
 *
 *----------------------------------------------------------------------------*/

void fisher_matrix(double *a, double *eval, double **evect, int ishuffle)
{
  int i,j,k,n, use_file = 0, i1;
  double dd = 0.01, *a_orig, x1, x2, x3, x4, x0;
  double **cov, **tmp, **tmp1, *atemp;
  int nrot, flag_normalize_eigenvalues = 0;
  long IDUM = -555;
  FILE *fp;
  char aa[1000];

  muh(0);
  exit(0);

  if(FISHER==3)
    flag_normalize_eigenvalues = 1;

  n = wpl.ncf;

  a_orig = dvector(1,n);
  cov = dmatrix(1,n,1,n);
  tmp = dmatrix(1,n,1,n);
  atemp = dvector(1,n);
  for(k=1;k<=n;++k) a_orig[k] = a[k];


  if(FISHER>1) 
    {
      use_file = 1;
      goto SKIP_COV_CALC2; // for reading in the evals/vects
      goto SKIP_COV_CALC;
    }

  for(i=1;i<=n;++i)
    a_orig[i] = a[i];

  wpl.reset_fred = 1;
  x0 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);
  printf("CHI %e\n",x0);

  for(i=1;i<=n;++i)
    for(j=i+1;j<=n;++j)
      {

	if(i==j) continue;
	// evaluate chi at i+,j+
	for(k=1;k<=n;++k) a[k] = a_orig[k];
	a[i] = a_orig[i]*(1+dd);
	a[j] = a_orig[j]*(1+dd);
	for(k=1;k<=n;++k) wpl.a[k] = a[k];
	wpl.reset_fred = 1;
	x1 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

	// evaluate chi at i+,j-
	for(k=1;k<=n;++k) a[k] = a_orig[k];
	a[i] = a_orig[i]*(1+dd);
	a[j] = a_orig[j]*(1-dd);
	for(k=1;k<=n;++k) wpl.a[k] = a[k];
	wpl.reset_fred = 1;
	x2 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

	// evaluate chi at i-,j+
	for(k=1;k<=n;++k) a[k] = a_orig[k];
	a[i] = a_orig[i]*(1-dd);
	a[j] = a_orig[j]*(1+dd);
	for(k=1;k<=n;++k) wpl.a[k] = a[k];
	wpl.reset_fred = 1;
	x3 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

	// evaluate chi at i-,j-
	for(k=1;k<=n;++k) a[k] = a_orig[k];
	a[i] = a_orig[i]*(1-dd);
	a[j] = a_orig[j]*(1-dd);
	for(k=1;k<=n;++k) wpl.a[k] = a[k];
	wpl.reset_fred = 1;
	x4 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

	cov[i][j] = 0.5*(x1-x2-x3+x4)/(4*dd*a_orig[i]*dd*a_orig[j]);
	tmp[i][j] = cov[i][j];

	cov[j][i] = cov[i][j];
	tmp[j][i] = cov[j][i];

	if(!use_file)
	  printf("CHI %d %d %e %e %e %e %e %e %e %e\n",
		 i,j,x1,x2,x3,x4,cov[i][j],dd*a_orig[i],a_orig[j]*dd,
		 (x1-x2-x3+x4));

      }

  // now for the daigonals
  for(i=1;i<=n;++i)
    {
      j = i;
      for(k=1;k<=n;++k) a[k] = a_orig[k];
      a[i] = a_orig[i]*(1+dd);
      for(k=1;k<=n;++k) wpl.a[k] = a[k];
      wpl.reset_fred = 1;
      x1 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

      for(k=1;k<=n;++k) a[k] = a_orig[k];
      a[i] = a_orig[i]*(1-dd);
      for(k=1;k<=n;++k) wpl.a[k] = a[k];
      wpl.reset_fred = 1;
      x2 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);

      cov[i][i] = 0.5*(x1-2*x0+x2)/(dd*dd*a_orig[i]*a_orig[i]);
      tmp[i][i] = cov[i][i];

      if(!use_file)
	printf("CHI %d %d %e %e %e\n",i,i,x1,x2,x0);
    }

 SKIP_COV_CALC:
  if(use_file)
    {
      fp = openfile(ARGV[4]);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  fscanf(fp,"%s %d %d %lf",aa,&i1,&i1,&cov[i][j]);
      fclose(fp);
      //temp
      fp = openfile("out2.diag");
      for(i=1;i<=-n;++i)
	fscanf(fp,"%s %d %d %lf",aa,&i1,&i1,&cov[i][i]);
      fclose(fp);
    }

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      {
	printf("COV %d %d %e\n",i,j,cov[i][j]);
      }
  fflush(stdout);

  tmp1 = dmatrix(1,n,1,1);

  gaussj(tmp,n,tmp1,1);
  jacobi(tmp,n,eval,evect,&nrot);

 SKIP_COV_CALC2:
  if(use_file)
    {
      fp = openfile(ARGV[4]);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  fscanf(fp,"%s %lf",aa,&evect[i][j]);
      fclose(fp);
      fp = openfile(ARGV[5]);
      for(i=1;i<=n;++i)
	fscanf(fp,"%s %lf",aa,&eval[i]);
      fclose(fp);
      for(i=1;i<=n;++i)
	eval[i] = fabs(eval[i]);
      
      if(flag_normalize_eigenvalues)
	{
	  normalize_eigenvalues(evect,eval,a);
	  for(i=1;i<=wpl.ncf;++i)
	    printf("EVALFAC %e\n",eval[i]);
	  exit(0);
	}

      if(ishuffle)
	{
	  for(i=1;i<=n;++i)
	    for(atemp[i]=0,j=1;j<=n;++j)
	      atemp[i] += sqrt(eval[j])*evect[j][i]*gasdev(&IDUM_MCMC);
	  for(i=1;i<=n;++i)
	    a[i] += atemp[i];
	}	  
      return ;

      // testing here
      //eval[1] *= 0.1;
      //eval[2] *= 0.1;
      //eval[12] *= 0.1;
      //eval[13] *= 0.1;
      fmuh(gasdev(&IDUM_MCMC));

      for(i=1;i<=wpl.ncf;++i)
	a_orig[i] = a[i];

      for(k=1;k<=wpl.ncf;++k)
	{
	  for(i=1;i<=wpl.ncf;++i)
	    a[i] = a_orig[i];

	  eval[k] *= 0.01;

	  for(i=1;i<=wpl.ncf;++i)
	    wpl.a[i] = a[i];
	  wpl.reset_inversion = 1;
	  wpl.reset_fred = 1;
	  i = check_lensing_parameters(a);
	  if(!i)
	    chi2_stellar_mass_function(a);
	  else
	    printf("CHISMF nan\n");
	  eval[k] *= 100;
	}

	  exit(0);
      return;
    }

  for(i=1;i<=n;++i)
    printf("EVAL %e\n",eval[i]);
  fflush(stdout);

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      printf("EVECT %e\n",evect[i][j]);
  fflush(stdout);

  for(i=1;i<=n;++i)
    eval[i] = fabs(eval[i]);
  muh(0);

  for(i=1;i<=n;++i)
    a[i] = a_orig[i];
  muh(1);

  //return;

  for(i=1;i<=n;++i)
    atemp[i] = sqrt(fabs(eval[i]));
  muh(2);
  
  for(i=1;i<=n;++i)
    for(a[i]=0,j=1;j<=n;++j)
      a[i] += atemp[j]*evect[j][i];
  muh(3);
  
  for(i=1;i<=n;++i)
    printf("PARAMp%d %e %e %e\n",i,a[i]+a_orig[i],a_orig[i],a[i]);

  for(i=1;i<=n;++i)
    atemp[i] = -sqrt(fabs(eval[i]));
  
  for(i=1;i<=n;++i)
    for(a[i]=0,j=1;j<=n;++j)
      a[i] += atemp[j]*evect[j][i];
  
  for(i=1;i<=n;++i)
    printf("PARAMm%d %e %e %e\n",i,a[i]+a_orig[i],a_orig[i],a[i]);

  exit(0);

}

/* This is to normalize the eigenvalues. Set a factor by which to multiple 
 * each eigenvalue such that when stepping 1-sigma, you get a delta chi^2=1
 */
void normalize_eigenvalues(double **evect, double *eval, double *a)
{
  double *fac;
  int i, j;

  fprintf(stdout,"norm> starting eigenvalue renormalization\n");
  fflush(stdout);

  // get the original chi2
  wpl.reset_fred = 1;
  chi2original_gl1 = chi2_stellar_mass_function(a) + 
    chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);
  chi2original_gl1 = chi2_stellar_mass_function(a) + 
    chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);
  printf("norm> chi2original = %e\n",chi2original_gl1);

  for(i=1;i<=wpl.ncf;++i)
    wpl.a[i] = a[i];

  eval_gl1 = dvector(1,wpl.ncf);
  evect_gl1 = dmatrix(1,wpl.ncf,1,wpl.ncf);

  for(i=1;i<=wpl.ncf;++i)
    eval_gl1[i] = eval[i];

  for(i=1;i<=wpl.ncf;++i)
    for(j=1;j<=wpl.ncf;++j)
      evect_gl1[i][j] = evect[i][j];

  fac = dvector(1,wpl.ncf);
  for(i=1;i<=wpl.ncf;++i)
    {
      eigeni_gl1 = i;
      fac[i] = zbrent(func_eigen_normalize,0.0,10.0,1.0E-2);
      printf("EIGNORM %d %e\n",i,fac[i]);
    }
  for(i=1;i<=wpl.ncf;++i)
    eval[i] = eval[i]*fac[i];

}

double func_eigen_normalize(double fac)
{
  static double *a;
  static int flag = 1;
  double atemp[100], x2, a_orig[100];
  int i, j, n;

  if(flag)
    {
      a = dvector(1,wpl.ncf);
      flag = 0;
    }

  n = wpl.ncf;

  for(i=1;i<=n;++i)
    a_orig[i] = wpl.a[i];

  for(i=1;i<=n;++i)
    atemp[i] = 0;
  atemp[eigeni_gl1] = sqrt(eval_gl1[eigeni_gl1])*fac;
  
  // get the perturbation to the position
  for(i=1;i<=n;++i)
    for(a[i]=0,j=1;j<=n;++j)
      a[i] += atemp[j]*evect_gl1[j][i];

  for(i=1;i<=-n;++i)
    if(fabs(a[i])>1.0E-10)printf("STEP %d %e %e %e\n",i,a[i],wpl.a[i],atemp[i]);

  // add to the original position
  for(i=1;i<=n;++i)
    a[i] += wpl.a[i];

  i = check_lensing_parameters(a);
  if(i) { 
    x2 = 1.0E+7;     
    printf("RANGE %d %e %e\n",i,a[abs(i)],a[abs(i)+HOD_PARAMS]);
    goto SKIPNORM1; }

  // reset the global params
  for(i=1;i<=n;++i)
    wpl.a[i] = a[i];

  // get the chi2
  wpl.reset_fred = 1;
  x2 = chi2_stellar_mass_function(a) + chi2_wp_lensing_wrapper(a) + chi2_lensing(a) + chi2_red_central_fraction(a);
  
  // put params back in place
  for(i=1;i<=n;++i)
    wpl.a[i] = a_orig[i];

 SKIPNORM1:
  printf("funcnorm> %d %e %e %e %e\n",eigeni_gl1,fabs(x2-chi2original_gl1) - 1,fac,x2,chi2original_gl1);
  return fabs(x2 - chi2original_gl1) - 1;
}

/*------------------------------------------------------------------------
 * POST-PROCESSING OF A FULL CHAIN 
 * ----------------------------------------------------------------------*/

double func_mean_hmass(double m)
{
  // want to return number of galaxies as a function of halo mass m
  m = exp(m);
  return N_cen(m)*dndM_interp(m)*m*m;
}

void chain_post_process()
{
  int i,j,k,n,i1,nm,ibuf=11;
  double *a, **red_sat_smf, **red_cen_smf, **blue_sat_smf, **blue_cen_smf,*xx,**red_mean_hmass,**blue_mean_hmass,
    **mhalo_blue, **mhalo_red, mass, *mass_array, ngalr, ngalb;
  FILE *fp, *fpout;
  char fname[1000], aa[1000];

  double mlo, mhi, dlogm;

  //set up parameter array
  a = dvector(1,wpl.ncf);
  // open the file for post processing
  fp = openfile(ARGV[2]);
  // number of elements
  n = filesize(fp);

  fprintf(stderr,"post_process> opening [%s]\n",ARGV[2]);
  fprintf(stderr,"post_process> processing [%d] elements\n",n);

  // set up the stellar mass array
  nm = 200;
  mlo = 1.0E+8;
  mhi = 1.0E+12;
  dlogm = log(mhi/mlo)/(nm-1);

  red_sat_smf = dmatrix(1,nm,1,n);
  red_cen_smf = dmatrix(1,nm,1,n);
  blue_sat_smf = dmatrix(1,nm,1,n);
  blue_cen_smf = dmatrix(1,nm,1,n);
  mhalo_blue = dmatrix(1,nm,1,n);
  mhalo_red =  dmatrix(1,nm,1,n);
  red_mean_hmass =  dmatrix(1,nm,1,n);
  blue_mean_hmass =  dmatrix(1,nm,1,n);

  // temp vector
  xx = dvector(1,n);

  // loop through all elements
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%s %d %d",aa,&i1,&i1);
      for(j=1;j<=wpl.ncf;++j)
	fscanf(fp,"%lf",&a[j]);
      fgets(aa,1000,fp);

      // set up globale
      for(j=1;j<=wpl.ncf;++j)
	wpl.a[j] = a[j];

      fprintf(stderr,"post_process> %d\n",i);

      // now loop through all masses
      for(j=1;j<=nm;++j)
	{
	  mass = exp((j-1)*dlogm)*mlo;
	  wpl.reset_fred = 1;
	  mhalo_blue[j][i] = ms_to_mhalo(mass,a);
	  mhalo_red[j][i] = ms_to_mhalo(mass,&(a[ibuf]));
	  
	  // get the stellar mass functions for red
	  BLUE_FLAG = 0;
	  wpl.reset_inversion = 1;
	  ms_to_mhalo_inversion(mhalo_red[j][i]);
	  set_up_hod_for_shmr(exp(-dlogm/2)*mass, exp(dlogm/2)*mass, a);
	  ngalr = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
	  red_cen_smf[j][i] = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
	  red_sat_smf[j][i] = ngalr-red_cen_smf[j][i];
	  red_mean_hmass[j][i] = qromo(func_mean_hmass,log(HOD.M_low),log(HOD.M_max),midpnt)/
	    red_cen_smf[j][i]/dlogm*log(10)*pow(HUBBLE,3.0);

	  // reset the tabulation of Ncen_xigm
	  set_up_hod_for_shmr(exp(-dlogm/2)*mass/2, exp(dlogm/2)*mass, a);
	  Nsat_xigm(mhalo_blue[j][i]);
	  
	  // get SMFs for blue gals
	  BLUE_FLAG = 1;
	  wpl.reset_inversion = 1;
	  set_up_hod_for_shmr(exp(-dlogm/2)*mass, exp(dlogm/2)*mass, a);
	  ngalb = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
	  blue_cen_smf[j][i] = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
      	  blue_sat_smf[j][i] = ngalb-blue_cen_smf[j][i];
	  blue_mean_hmass[j][i] = qromo(func_mean_hmass,log(HOD.M_low),log(HOD.M_max),midpnt)/
	    blue_cen_smf[j][i]/dlogm*log(10)*pow(HUBBLE,3.0);

	}
    }

  // now that arrays are in place, find 68% values for each desired quantity
  sprintf(fname,"stats_z%d.fred_halomass",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,blue_cen_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,blue_cen_smf[i][n/2],
	     blue_cen_smf[i][(int)(n*0.16)],blue_cen_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.meanhalomass",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,blue_mean_hmass[i]);
      sort(n,red_mean_hmass[i]);
      fprintf(fpout,"%e %e %e %e ",exp((i-1)*dlogm)*mlo,blue_mean_hmass[i][n/2],
	     blue_mean_hmass[i][(int)(n*0.16)],blue_mean_hmass[i][(int)(n*0.84)]);
      fprintf(fpout,"%e %e %e\n",red_mean_hmass[i][n/2],
	     red_mean_hmass[i][(int)(n*0.16)],red_mean_hmass[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.blue_cen_smf",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,blue_cen_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,blue_cen_smf[i][n/2],
	     blue_cen_smf[i][(int)(n*0.16)],blue_cen_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.blue_sat_smf",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,blue_sat_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,blue_sat_smf[i][n/2],
	     blue_sat_smf[i][(int)(n*0.16)],blue_sat_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.red_cen_smf",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,red_cen_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,red_cen_smf[i][n/2],
	     red_cen_smf[i][(int)(n*0.16)],red_cen_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.red_sat_smf",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,red_sat_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,red_sat_smf[i][n/2],
	     red_sat_smf[i][(int)(n*0.16)],red_sat_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.fsat_red",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = red_sat_smf[i][j]/(red_sat_smf[i][j]+red_cen_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.fsat_blue",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = blue_sat_smf[i][j]/(blue_sat_smf[i][j]+blue_cen_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.fq_sat",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = red_sat_smf[i][j]/(blue_sat_smf[i][j]+red_sat_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.fq_cen",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = red_cen_smf[i][j]/(blue_cen_smf[i][j]+red_cen_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);


  sprintf(fname,"stats_z%d.fq_all",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = (red_cen_smf[i][j]+red_sat_smf[i][j])/
		 (blue_cen_smf[i][j]+red_cen_smf[i][j]+blue_sat_smf[i][j]+red_sat_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.fsat_all",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = (red_sat_smf[i][j]+blue_sat_smf[i][j])/
		 (blue_cen_smf[i][j]+red_cen_smf[i][j]+blue_sat_smf[i][j]+red_sat_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);


  sprintf(fname,"stats_z%d.shmr_blue",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = mhalo_blue[i][j];
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"stats_z%d.shmr_red",wpl.iz);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = mhalo_red[i][j];
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);
  exit(0);
}

/*------------------------------------------------------------------------
 * POST-PROCESSING OF CHAINS TO GET ERROR IN RED/BLUE LETTER
 * ----------------------------------------------------------------------*/

void chain_post_process_central_red_error()
{
  int i,j,k,n,i1,nm,ibuf=11,i2;
  double *a, **red_sat_smf, **red_cen_smf, **blue_sat_smf, **blue_cen_smf,*xx,
    **mhalo_blue, **mhalo_red, *mass_array, ngalr, ngalb, *temp, mstar;
  FILE *fp, *fpout;
  char fname[1000], aa[1000];

  double mlo, mhi, dlogm;

  int nhalo;
  double *mass;
  long IDUM=-555;

  // get group galaxies
  fp = openfile("central_galaxies_red.dat");
  nhalo = filesize(fp);
  mass = dvector(1,nhalo);
  temp = dvector(1,nhalo);

  for(i=1;i<=nhalo;++i) {
    fscanf(fp,"%lf",&mass[i]); fgets(aa,1000,fp); }
  fclose(fp);
  for(i=1;i<=nhalo;++i)
    mass[i] = pow(10.0,mass[i])*HUBBLE*1.1; //conver to h-inverse units (and 200c to 220m)




  //set up parameter array
  a = dvector(1,wpl.ncf);
  // open the file for post processing
  fp = openfile(ARGV[2]);
  // number of elements
  n = filesize(fp);


  fprintf(stderr,"post_process> opening [%s]\n",ARGV[2]);
  fprintf(stderr,"post_process> processing [%d] elements\n",n);

  BLUE_FLAG = 0;

  // loop through all elements
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%s %d %d",aa,&i1,&i2);
      for(j=1;j<=wpl.ncf;++j)
	fscanf(fp,"%lf",&a[j]);
      fgets(aa,1000,fp);

      if(i1<1000)continue;

      // set up global
      for(j=1;j<=wpl.ncf;++j)
	wpl.a[j] = a[j];

      // for each halo in catalog, randomly sample central red galaxy
      for(j=1;j<=nhalo;++j)
	{
	  mstar = log10(ms_to_mhalo_inversion(mass[j]));
	  temp[j] = gasdev(&IDUM)*wpl.a[6+(1-BLUE_FLAG)*11]+(mstar);
	}
      sort(nhalo,temp);
      printf("MEDIAN %d %f\n",i,temp[nhalo/2]);
    }
  exit(0);
}

/*------------------------------------------------------------------------
 * POST-PROCESSING OF CHAINS OR ELEMENTS IN THE CHAIN
 * ----------------------------------------------------------------------*/

void chain_post_process_single(double *a)
{
  int i,j,k,n, i1, i2;
  //  double *a;
  char aa[1000];
  FILE *fp, *fpout;
  

  double mass, mlo, mhi, dlogm, fsatr, fsatb, fred, ngal, ngalb, ngalr, mhalor, mhalob, nsatr, nsatb, ncenr, ncenb;
  int nm=200, ibuf;

  n = wpl.ncf;
  ibuf = HOD_PARAMS;
  sprintf(aa,"stats1.z%d",wpl.iz);
  fpout = fopen(aa,"w");

  muh(ibuf);
  /*
  wpl.ncf = n = 27;
  a = dvector(1,n);
  HOD_PARAMS = ibuf = 11;

  fp = openfile(ARGV[2]);
  fscanf(fp,"%s %d %d",aa,&i1,&i2);
  for(i=1;i<=n;++i)
    fscanf(fp,"%lf",&a[i]);
  fclose(fp);

  printf("POST PROCESS FILE: [%s]\n",ARGV[2]);
  for(i=1;i<=n;++i)
    printf("a[%02d]= %e\n",i,a[i]);
  for(i=1;i<=n;++i)
    wpl.a[i] = a[i];
  */

  mlo = 1.0E+8;
  mhi = 1.0E+12;
  dlogm = log(mhi/mlo)/(nm-1);

  for(i=1;i<=nm;++i)
    {
      mass = exp((i-1)*dlogm)*mlo;
      mhalob = ms_to_mhalo(mass,a);
      mhalor = ms_to_mhalo(mass,&(a[ibuf]));

      
      BLUE_FLAG = 0;
      wpl.reset_inversion = 1;
      ms_to_mhalo_inversion(mhalor);
      set_up_hod_for_shmr(exp(-dlogm/2)*mass, exp(dlogm/2)*mass, a);
      ngalr = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      nsatr = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      ncenr = ngalr - nsatr;
      fsatr = nsatr/ngalr;

      // reset the tabulation of Ncen_xigm
      set_up_hod_for_shmr(exp(-dlogm/2)*mass/2, exp(dlogm/2)*mass, a);
      Nsat_xigm(mhalob);

      BLUE_FLAG = 1;
      wpl.reset_inversion = 1;
      set_up_hod_for_shmr(exp(-dlogm/2)*mass, exp(dlogm/2)*mass, a);
      ngalb = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      nsatb = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      ncenb = ngalb - nsatb;
      fsatb = nsatb/ngalb;
      

      //printf("HODRED %e %e %e %e %e %e %e\n",mass,mhalor,HOD.M_low,
      //     HOD.M_min,red_central_fraction(HOD.M_min,a), N_cen(HOD.M_min),ncenr);

      fprintf(fpout,"STATS1 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
	     mass, mhalob, mhalor, fsatb, fsatr, nsatb/(ngalb+ngalr), nsatr/(ngalb+ngalr), 
	     ngalb/dlogm*log(10)*pow(HUBBLE,3.0), ngalr/dlogm*log(10)*pow(HUBBLE,3.0),
	     ncenr/(ncenr+ncenb),nsatr/(nsatr+nsatb),
	     nsatb/dlogm*log(10)*pow(HUBBLE,3.0), nsatr/dlogm*log(10)*pow(HUBBLE,3.0),
	     ncenb/dlogm*log(10)*pow(HUBBLE,3.0), ncenr/dlogm*log(10)*pow(HUBBLE,3.0));
    }
  fclose(fpout);

  mlo = 1.0E+10;
  mhi = 1.0E+14;
  dlogm = log(mhi/mlo)/(nm-1);

  sprintf(aa,"stats2.z%d",wpl.iz);
  fpout = fopen(aa,"w");

  for(i=1;i<=nm;++i)
    {
      mass = exp((i-1)*dlogm)*mlo;
      fred = red_central_fraction(mass, a);
      fprintf(fpout,"STATS2 %e %e\n",mass,fred);
    }

}

/*******************************************************************
 * -----------------------------------------------------------------
 * do CLF post-processing
 * -----------------------------------------------------------------
 *******************************************************************/

void chain_post_process_clf(double *a)
{
  int i,j,k,n, i1, i2;
  //  double *a;
  char aa[1000];
  FILE *fp, *fpout;
  

  double mass, mlo, mhi, dlogm, fsatr, fsatb, fred, ngal, ngalb, ngalr, mhalor, mhalob, nsatr, nsatb, ncenr, ncenb, dm;
  int nm=20, ibuf;
  double ngalred,ngalblue,ngaltot,mgal;

  n = 100;

  for(j=120;j<=150;++j)
    {
      mass = pow(10.0,j/10.0);

      mlo = 10.0;
      mhi = 10.2;
      dm = (mhi-mlo)/(n-1);

      // get the red satellites galaxies
      BLUE_FLAG = 0;
      wpl.reset_inversion = 1;
      ms_to_mhalo_inversion(mass);

      ngaltot = 0;
      for(i=1;i<=n;++i)
	{
	  mgal = (i-1)*dm + mlo;
	  set_up_hod_for_shmr(pow(10.0,mgal-dm/2),pow(10.0,mgal+dm/2),wpl.a);
	  ngal = Nsat_xigm(mass);
	  ngaltot += ngal;
	}
      ngalred = ngaltot;

      // get the red satellites galaxies
      BLUE_FLAG = 1;
      wpl.reset_inversion = 1;
      ms_to_mhalo_inversion(mass);

      ngaltot = 0;
      for(i=1;i<=n;++i)
	{
	  mgal = (i-1)*dm + mlo;
	  set_up_hod_for_shmr(pow(10.0,mgal-dm/2),pow(10.0,mgal+dm/2),wpl.a);
	  ngal = Nsat_xigm(mass);
	  ngaltot += ngal;
	}
      ngalblue = ngaltot;
      printf("REDFRAC %e %e %e %e\n",mass,ngalred/(ngalred+ngalblue),ngalred,ngalblue);
    }
  exit(0);
	     
}

/*******************************************************************
 * -----------------------------------------------------------------
 * do redshift-dependent post-processing
 * -----------------------------------------------------------------
 *******************************************************************/

void quenching_halo_mass()
{

}


/*******************************************************************
 * -----------------------------------------------------------------
 * do a chi^2 minimization
 * -----------------------------------------------------------------
 *******************************************************************/

void chi2_lensing_minimization()
{
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min;

  n = wpl.ncf;

  a = dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  initial_lensing_values(a,pp,yy);

  if(POWELL) 
    {
      if(OUTPUT)printf("lensing_min> starting powell.\n");
      powell(a,pp,n,FTOL,&niter,&chi2min,chi2_lensing_total);
      chi2min = chi2_lensing_total(a);
    }
  else
    {
      if(OUTPUT)printf("lensing_min> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_lensing_total,&niter);
      for(i=1;i<=n;++i)a[i]=pp[1][i];
      chi2min = chi2_lensing_total(a);
    }	

  printf("POWELL %e ",chi2min);
  for(i=1;i<=n;++i)printf("%e ",a[i]);
  printf(" \n");
  exit(0);
}

double chi2_lensing_total(double *a)
{
  int i;
  double chi2, chi2smf, chi2wp, chi2fred, chi2lensing;
  static int niter=0;

  wpl.reset_inversion = 1;
  wpl.reset_fred = 1;

  for(i=1;i<=wpl.ncf;++i)
    wpl.a[i] = a[i];

  if(check_lensing_parameters(a))
    {
      chi2 = 1.0E+7;
      goto SKIPCHI;
    }

  chi2smf = chi2_stellar_mass_function(a);
  chi2wp = chi2_wp_lensing_wrapper(a);
  printf("> %d %d %d\n",HOD.pdfc,HOD.pdfs,EXCLUSION);
  chi2lensing = chi2_lensing(a);
  printf("> %d %d %d\n",HOD.pdfc,HOD.pdfs,EXCLUSION);
  chi2fred = chi2_red_central_fraction(a);
  chi2 = chi2wp + chi2lensing + chi2smf + chi2fred;

 SKIPCHI:

  printf("ITER %d %e ",++niter,chi2);
  for(i=1;i<=wpl.ncf;++i)
    printf("%e ",wpl.a[i]);
  printf("\n");
  fflush(stdout);
  return chi2;
}

void initial_lensing_values(double *a, double **pp, double *yy)
{
  int n, i, j;
  FILE *fp;
  char aa[1000];
  double d[1000];

  n = wpl.ncf;
  muh(wpl.ncf);

  fp = openfile(ARGV[2]);
  fscanf(fp,"%s %d %d",aa,&i,&j);
  for(i=1;i<=wpl.ncf;++i)
    fscanf(fp,"%lf",&a[i]);
  fclose(fp);
  
  for(i=1;i<=n;++i)
    d[i] = 0.1*a[i];

  d[1] = 0.01*a[1];
  d[2] = 0.01*a[2];
  d[12] = 0.01*a[12];
  d[13] = 0.01*a[13];

  /*
  for(i=1;i<=n;++i)
    a[i] = a[i] + (1-2*drand48())*d[i];
  */

  for(i=1;i<=n;++i)
    printf("a[%02d]= %e %e\n",i,a[i],d[i]);

  if(POWELL)
    {
      for(i=1;i<=wpl.ncf;++i)
	{
	  for(j=1;j<=wpl.ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j){ pp[i][j]+=d[j]; }
	    }
	}
    }
  else
    {
      for(j=1;j<=wpl.ncf;++j)
	pp[1][j]=a[j];
      yy[1]=chi2_lensing_total(a);
    
      for(i=1;i<=wpl.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	    yy[i+1]=chi2_lensing_total(a);	  
	  for(j=1;j<=wpl.ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[wpl.ncf]-=d[wpl.ncf];
    }

}
