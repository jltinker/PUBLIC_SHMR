#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* External functions from wp_minimization.c
 */
void wp_input(void);
double mcmc_initialize(double *a, double **cov1, double *avg1, double *start_dev);
double chi2rsd(double *a);
double chi2rsdq(double *a);


double chi2rsd(double *a)
{
  return 1;
}
double chi2rsdq(double *a)
{
  return 1;
}

/* Internal functions.
 */
double chi2_wp_wrapper(double *a);
void mcmc_restart2(double *start_dev, int np);
int mcmc_restart3(double **chain, int n, double *chi2_prev, int *iweight);
void qp_input_mcmc();
void output_mcmc_quants(int);

int USE_IWEIGHT = 0;
double omz0, s8z0, sigv0;

/**************************************
 * 
 * RESART OPTIONS
 * 
 * 1 - read in an continue
 * 2 - use exponential decay on chi2
 * 4 - do not update covariance matrix
 *
 */

/******************************************************************
 *
 * HOD.free[] also controls which variables will be held constant/vary
 * during MCMC minimization. Since this routine will also so z-space
 * minimization if requested, indices>6 are cosmological.
 *
 *  i     variable
 * ---    --------
 * [1] ->  M_min
 * [2] ->  M1
 * [3] ->  alpha
 * [4] ->  M_cut
 * [5] ->  sigmaM
 * [6] ->  CVIR_FAC
 * [7] ->  MaxCen (or M_cen_max)
 * [8] ->  M_sat_break
 * [9] ->  alpha1
 *
 * [10]->  OMEGA_M
 * [11]->  SIGMA_8
 * [12]->  VBIAS
 * [13]->  VBIAS_C
 * [14]->  GAMMA
 * [15]->  SPECTRAL_INDX
 *
 * [0] -> The galaxy_density will be considered data with errors on it,
 *         and therefore no variable will be fixed by the galaxy density.
 * 
 */
void mcmc_minimization()
{
  FILE *fpmc;
  char fname[1000];
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,xbias,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=10000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b;


  N_HOD_PARAMS = 9;
  sigv0 = SIGV;
  omz0 = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  s8z0 = SIGMA_8;

  opar=dvector(1,100);

  MCMC=Task.MCMC;
  if(ARGC>2)
    IDUM = -1*atoi(ARGV[2]);
  printf("IDUM= %d %d\n",(int)IDUM, MCMC);

  pcheck=calloc(ptot,sizeof(int));

  wp_input();
  if(MCMC>2)qp_input_mcmc();
  
  Work.imodel=2;
  Work.chi2=1;

  sprintf(fname,"%s.MCMC_%d", Task.root_filename, abs(IDUM));
  fprintf(stderr,"mcmc> output filename= [%s]\n",fname);
  fpmc = fopen(fname,"w");

  /*
  OUTPUT=0;
  */

  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  for(n=0,i=1;i<100;++i)
    {
      n+=HOD.free[i];
      /* if(i>N_HOD_PARAMS && HOD.free[i])MCMC=3;*/
      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  wp.ncf=n;

  if(HOD.free[0])
    {
      wp.ngal = GALAXY_DENSITY;
      wp.ngal_err = 0.1*wp.ngal;
      FIX_PARAM = 0;
    }

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

  if(RESTART)
    {
      niter = mcmc_restart3(chain,n,&chi2prev,iweight);
      if(niter < NSTEP)
	{
	  if(ThisTask==0)
	    fprintf(stderr,"Not enough points in restart chain: %d<=%d\n",niter,NSTEP);
	  exit(0);
	}
      for(i=1;i<=n;++i)
	aprev[i] = chain[niter][i];
      goto RESTART_POINT;
    }

  chi2prev=mcmc_initialize(a,cov1,avg1,start_dev);
  niter++;
  for(i=1;i<=n;++i)
    {
      aprev[i] = a[i];
      chain[1][i] = a[i];
    }

  pcnt=0;
  pcheck[pcnt]=1;

  if(RESTART==2)
    {
      mcmc_restart2(start_dev,n);
    }

  stepfac=wpl.stepfac_burn;
  if(wpl.stepfac_burn>0)stepfac = wpl.stepfac_burn;
  while(niter<NSTEP)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  //stepfac = stepfac*pow(0.9,5-j);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8Z0 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  if(HOD.free[++i])HUBBLE   = a[++j];
	  if(HOD.free[++i])OMEGA_B   = a[++j];
	  SIGMA_8 = SIGMA_8Z0*growthfactor(REDSHIFT);
	  fmuh(SIGMA_8);
	  dndM_interp(1.0E+13);
	}
      if(VBIAS_C<0)continue;

      /* Draw random value of cvir from prior.
       */
      /* if(CVIR_FAC<0.3 || CVIR_FAC>1.2)continue; */
      /* CVIR_FAC = 0.9*drand48()+0.3;  */
      /* GAMMA = gasdev(&IDUM)*0.02 + 0.15; */

      chi2=chi2_wp_wrapper(a);

      if(!ThisTask){
	printf("TRY %d ",++count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("%e\n",chi2);fflush(fpmc);
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

      // output the stats
      //output_mcmc_quants(niter);
      if(!ThisTask){
	xbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	fprintf(fpmc, "%d %d ",niter,count);
	for(i=1;i<=n;++i)
	  fprintf(fpmc,"%e ",a[i]);
	fprintf(fpmc,"%e %e %e %e\n",chi2, HOD.M_min, fsat, xbias);fflush(fpmc);
	printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
	       number_weighted_central_mass(),
	       qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);

	fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	printf("FSAT %d %e %e %e %e\n",niter,fsat,HOD.M_min,HOD.sigma_logM,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
      }

    }
 RESTART_POINT:

  stepfac=1.6/sqrt(n);
  pcnt=-1;
  t0 = second();

  NSTEP = niter;

  while(niter<imax_chain)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  //stepfac=1.6/sqrt(n);
	  //stepfac=0.6/sqrt(n);
	  //	  stepfac = stepfac*pow(0.9,6-j);
	  stepfac = 0.25;
	  stepfac = 0.5;
	  stepfac=1.6/sqrt(n);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac=1.6/sqrt(n);
      if(wpl.stepfac>=0)stepfac = wpl.stepfac;
      //stepfac = 0;

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

      /* We seem to be having a problem with this.
       * So, broadcast the model params from the root processor.
       */
#ifdef PARALLEL      
      MPI_Bcast(&a[1],n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
#endif

      // CHeck that the chi2 for the last point in the restart chain
      // is recovered.
      /*
      if(RESTART && !count)
	for(i=1;i<=n;++i)
	  a[i] = aprev[i];
      */
      if(RESTART==5)
	{
	  j = count+1;
	  if(j>4000)exit(0);
	  for(i=1;i<=n;++i)
	    a[i] = chain[j][i];
	}

      /*
      if(!ThisTask)
	for(i=1;i<=n;++i)
	  {
	    printf("COV %d %d %e ",count,i,sqrt(eval[i]));
	    for(j=1;j<=n;++j)
	      printf("%e ",evect[j][i]);
	    printf("\n");
	  }
      */

      /* Using only variances
       */
      //for(i=1;i<=n;++i) 
      //	a[i] = aprev[i] + gasdev(&IDUM)*sqrt(tmp[i][i])*stepfac;

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8Z0 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  if(HOD.free[++i])HUBBLE   = a[++j];
	  if(HOD.free[++i])OMEGA_B   = a[++j];
	  SIGMA_8 = SIGMA_8Z0*growthfactor(REDSHIFT);
	  dndM_interp(1.0E12);
	}
      if(VBIAS_C<0)continue;

      chi2 =chi2_wp_wrapper(a);

      if(RESTART==2)
	chi2*=(1+exp(-count/100.0));
      if(RESTART==3)
	chi2*=(1+exp(-count/20.0));

      tprev = t0;
      t0 = second();
      ++count;
      if(!ThisTask) {
	printf("TRY %d ",count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	if(RESTART==2) {
	  printf("%e %e %.2f\n",chi2,chi2/(1+exp(-count/100.0)),
		 timediff(tprev,t0));fflush(stdout); }
	else {
	  printf("%e %.2f\n",chi2,
		 timediff(tprev,t0));fflush(stdout); }
      }
      if(0) {
	printf("CPU%02d %d ",ThisTask,count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	if(RESTART==2) {
	  printf("%e %e %.2f\n",chi2,chi2/(1+exp(-count/100.0)),
		 timediff(tprev,t0));fflush(stdout); }
	else {
	  printf("%e %.2f\n",chi2,
		 timediff(tprev,t0));fflush(stdout); }
      }

      pcheck[pcnt]=0;
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
      pcheck[pcnt]=1;

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

      //output_mcmc_quants(niter);

      if(!ThisTask){
	xbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	fprintf(fpmc, "%d %d ",niter,count);
	for(i=1;i<=n;++i)
	  fprintf(fpmc,"%e ",a[i]);
	fprintf(fpmc,"%e %e %e %e\n",chi2, HOD.M_min, fsat, xbias);fflush(fpmc);
	printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
	       number_weighted_central_mass(),
	       qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);

	fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	printf("FSAT %d %e %e %e %e\n",niter,fsat,HOD.M_min,HOD.sigma_logM,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
      }

    }
}

void output_mcmc_quants(int niter)
{
  double x1,x2,x3,x4,x5,dr,r;
  int j;
  FILE *fp1;
  char fname[100];

  sprintf(fname,"wp.%d",niter);
  fp1=fopen(fname,"w");
  dr=(log(70.0)-log(0.01))/49.0;
  for(j=0;j<50;++j)
    {
      r=exp(j*dr+log(0.01));
      x1=one_halo_real_space(r);
      x2=two_halo_real_space(r);
      x3=projected_xi(r);
      x4 = projected_xi_1halo(r);
      x5 = projected_xi_rspace(r);
      fprintf(fp1,"%f %e %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4,x5);
      fflush(fp1);
    }
  fclose(fp1);

  KAISER=1;
  RESET_ZSPACE++;
  RESET_KAISER++;
  x1 = qromo(func_galaxy_bias, log(HOD.M_low), log(HOD.M_max), midpnt)/GALAXY_DENSITY;
  BETA = pow(OMEGA_Z,GAMMA)/x1;
  fprintf(stdout,"BETA= %f %f %f\n",BETA, OMEGA_Z, x1);
  fflush(stdout);
  xi_multipoles(r,&x1,&x2);
}

double chi2_wp_wrapper(double *a)
{
  static int flag=1, imcut =0;
  static double *b, chi2, chi2b, gbias;
  int i,j;
  double mcut, m1;

  m1 = HOD.M1;

  if(flag)
    {
      b=dvector(1,100);
      flag=0;
    }

  for(j=0,i=1;i<=N_HOD_PARAMS;++i) {
    if(HOD.free[i]) { 
      if(a[++j]<=0) { printf("NEG %d %d %e\n",i,j,a[j]);  return(1.0E7); } }
    if(HOD.free[i]) {
      ++j; }
  }

  i=0;j=0;
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_min */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]); m1 = b[j]; } /* M1 */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]); imcut = 1; mcut = b[j]; } /* M_cut */
  if(HOD.free[++i]){j++;b[j]=a[j];} /* sigma_logM */
  if(HOD.free[++i]){j++;b[j]=CVIR_FAC=a[j];}           /* cvir_fac */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* MaxCen */
  //if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_sat_break */
  //if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha1 */

  if(imcut && mcut/m1 < 0.0001)return 1.0E+7;
  
  // get the chi2 for wp
  chi2 = chi2_wp(b);

  // we've already set the cosmo parameters, now we've set the HOD params.
  gbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  OMEGA_Z = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  BETA = pow(OMEGA_Z,GAMMA)/gbias;
  SIGV = sigv0*pow(OMEGA_Z/omz0,GAMMA)*(SIGMA_8Z0/s8z0);
  printf("BETA %f bias= %f\n",BETA,gbias);

  if(MCMC>2) // get the RSD chi2
    {
      if(MCMC==2)
	chi2b = chi2rsd(b);
      if(MCMC==3)
	chi2b = chi2rsdq(b);
    }
  return(chi2+chi2b);
}

double mcmc_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0,k;
  double x1,x2,omega_m, deltax=1.0E-3,a0[100],fisher;
  long IDUM = -556;
  FILE *fp;

  omega_m = 1;

  i=0;j=0;
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_min/omega_m);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M1/omega_m);start_dev[j]=0.001; } //.0005
  if(HOD.free[++i]){ a[++j]=HOD.alpha;start_dev[j]=0.03; } //.005
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_cut/omega_m);start_dev[j]=0.01; } //.001
  if(HOD.free[++i]){ a[++j]=HOD.sigma_logM;start_dev[j]=0.01; }
  if(HOD.free[++i]){ a[++j]=CVIR_FAC;start_dev[j]=0.02; }
  if(HOD.pdfc==7) {
    if(HOD.free[++i])a[++j]=log10(HOD.M_cen_max/omega_m); start_dev[j]=0.001; }
  else {
    if(HOD.free[++i])a[++j]=HOD.MaxCen; start_dev[j]=0.02; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_sat_break/omega_m);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=HOD.alpha1;start_dev[j]=0.02; }

  if(MCMC>1)
    {
      if(HOD.free[++i])a[++j]=OMEGA_M;
      if(HOD.free[++i])a[++j]=SIGMA_8Z0;
      if(HOD.free[++i])a[++j]=VBIAS;
      if(HOD.free[++i])a[++j]=VBIAS_C;
      if(HOD.free[++i])a[++j]=GAMMA;
      if(HOD.free[++i])a[++j]=SPECTRAL_INDX;
      if(HOD.free[++i])a[++j]=HUBBLE;
      if(HOD.free[++i])a[++j]=OMEGA_B;
    }
  
  if(ARGC>3)
    {
      fp = openfile(ARGV[3]);
      fscanf(fp,"%d %d",&i,&j);
      for(i=1;i<=wp.ncf;++i)
	fscanf(fp,"%lf",&a[i]);
      fclose(fp);
    }
  

  if(!ThisTask)
    {
      printf("INITIAL VALUES: ");
      for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
      printf("\n");
    }

  for(i=1;i<=wp.ncf;++i)
    {
      avg1[i]=a[i];
      for(j=1;j<=wp.ncf;++j)
	cov1[i][j]=a[i]*a[j];
    }

  if(MCMC>1)
    {
      RESET_COSMOLOGY++;
      j=0;
      for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
      i=N_HOD_PARAMS;
      if(HOD.free[++i]){ OMEGA_M = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ SIGMA_8Z0 = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS   = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS_C = a[++j]; start_dev[j] = 0.02; } 
      if(HOD.free[++i]){ GAMMA   = a[++j]; start_dev[j] = 0.015; } 
      if(HOD.free[++i]){ SPECTRAL_INDX    = a[++j]; start_dev[j] = 0.02; }
      if(HOD.free[++i]){ HUBBLE    = a[++j]; start_dev[j] = 0.02; }
      if(HOD.free[++i]){ OMEGA_B    = a[++j]; start_dev[j] = 0.02; }
      SIGMA_8 = SIGMA_8Z0*growthfactor(REDSHIFT);
    }

  x1=chi2_wp_wrapper(a);
  
  // get the start dev values from the diagonal elements of the Fisher matrix
  // get the ii chi2 values
  for(i=1;i<=wp.ncf;++i)
    a0[i] = a[i];
  for(k=1;k<=-wp.ncf;++k)
    {
      a[k] = a0[k]*(1+deltax);
      RESET_COSMOLOGY++;
      j=0;
      for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
      i=N_HOD_PARAMS;
      if(HOD.free[++i]){ OMEGA_M = a[++j]; } 
      if(HOD.free[++i]){ SIGMA_8Z0 = a[++j]; } 
      if(HOD.free[++i]){ VBIAS   = a[++j];} 
      if(HOD.free[++i]){ VBIAS_C = a[++j]; } 
      if(HOD.free[++i]){ GAMMA   = a[++j];  } 
      if(HOD.free[++i]){ SPECTRAL_INDX    = a[++j];  }
      if(HOD.free[++i]){ HUBBLE    = a[++j];  }
      if(HOD.free[++i]){ OMEGA_B    = a[++j]; }

      x1 = chi2_wp_wrapper(a);
      fisher = 2*x1/(deltax*deltax*a0[k]*a0[k]);
      start_dev[k] = sqrt(1/fisher)/a0[k];
      printf("FISHER: %d %e %e\n",k,start_dev[k],x1);
    }

  // output_mcmc_quants(1);

  x2=100; // set high to get ball rolling.

  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e\n",x1+x2);fflush(stdout);
    printf("INITIAL CHI2: %e %e\n",x1,x2);
    fflush(stdout);
  }
  return(x1+x2);
}

/* This is to look at a chain and get the variances in each parameter.
 */
void mcmc_restart2(double *start_dev, int np)
{
  int n,i,j,k,i1,i2;
  FILE *fp;
  char aa[100];
  float xbar[10],xsqr[10],x;

  fp = openfile(RESTART_FILE);
  n = filesize(fp);

  for(i=0;i<np;++i)
    xbar[i] = xsqr[i] = 0;

  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%s %d %d",aa,&i1,&i2);
      for(j=0;j<np;++j)
	{
	  fscanf(fp,"%f",&x);
	  xbar[j] += x;
	  xsqr[j] += x*x;
	}
      fgets(aa,100,fp);
    }
  for(i=0;i<np;++i)
    {
      xbar[i]/=n;
      xsqr[i]/=n;
      xsqr[i] = sqrt(xsqr[i] - xbar[i]*xbar[i]);
      start_dev[i+1] = xsqr[i];
      if(!ThisTask)
	fprintf(stdout,"RESTART_DEV %f %f\n",xbar[i],xsqr[i]);
    }
}

int mcmc_restart3(double **chain, int n, double *chi2_prev, int *iweight)
{
  FILE *fp;
  char aa[100];
  int niter,i,j,i1,i2,iprev;
  double x,*a,chi2;

  fp = openfile(RESTART_FILE);
  niter = filesize(fp);

  a = dvector(1,n);

  fscanf(fp,"%s %d %d",aa,&i1,&i2);
  rewind(fp);
  iprev = i2 - 1;
  for(i=1;i<=niter;++i)
    {
      fscanf(fp,"%s %d %d",aa,&i1,&i2);
      if(USE_IWEIGHT)
	iweight[i] = i2 - iprev;
      else
	iweight[i] = 1;
      iprev =  i2;
      for(j=1;j<=n;++j)
	fscanf(fp,"%lf",&chain[i][j]);
      fscanf(fp,"%lf",&x);
      fgets(aa,100,fp);
    }
  if(RESTART==2 || RESTART==3)
    *chi2_prev = 20000*x; // set it to automatically take the first element
  else
    *chi2_prev = x;

  /* Normalize all the masses by OMEGA_M
   */
  for(i=1;i<=-niter;++i)
    {
      chain[i][1] -= log10(chain[i][4]);
      chain[i][3] -= log10(chain[i][4]);
    }
  return niter;
}

void qp_input_mcmc()
{
  char fname[1000],aa[1000];
  FILE *fp;
  int i,j,k;
  float x1,x2,x3;
  double **tmp, **tmp2;

  sprintf(fname,"monopole_mockmean.dat");
  fp = openfile(fname);
  rsdm.np = filesize(fp);
  rsdm.covar = dmatrix(1,rsdm.np,1,rsdm.np);
  rsdm.r = dvector(1,rsdm.np);
  rsdm.x = dvector(1,rsdm.np);
  rsdm.e = dvector(1,rsdm.np);

  for(j=1;j<=rsdm.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      rsdm.x[j] = x2;
      rsdm.r[j] = x1;
      rsdm.e[j] = x3;
    }
  fclose(fp);


  sprintf(fname,"monopole_mockmean.covar");
  fp = openfile(fname);
  for(j=1;j<=rsdm.np;++j)
    for(k=1;k<=rsdm.np;++k)
      fscanf(fp,"%lf",&rsdm.covar[j][k]);
  fclose(fp);

  // inver teh matrix
  tmp=dmatrix(1,rsdm.np,1,1);
  tmp2=dmatrix(1,rsdm.np,1,rsdm.np);
  for(i=1;i<=rsdm.np;++i)
    for(j=1;j<=rsdm.np;++j)
      tmp2[i][j]=rsdm.covar[i][j];
  gaussj(tmp2,rsdm.np,tmp,1);
  for(i=1;i<=rsdm.np;++i)
    for(j=1;j<=rsdm.np;++j)
      rsdm.covar[i][j]=tmp2[i][j];
  free_dmatrix(tmp,1,rsdm.np,1,1);
  free_dmatrix(tmp2,1,rsdm.np,1,rsdm.np);

  /*-- Now do the quadrupole
   */

  sprintf(fname,"Q_mockmean.dat");
  fp = openfile(fname);
  rsdq.np = filesize(fp);
  rsdq.covar = dmatrix(1,rsdq.np,1,rsdq.np);
  rsdq.r = dvector(1,rsdq.np);
  rsdq.x = dvector(1,rsdq.np);
  rsdq.e = dvector(1,rsdq.np);

  for(j=1;j<=rsdq.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      rsdq.x[j] = x2;
      rsdq.r[j] = x1;
      rsdq.e[j] = x3;
    }
  fclose(fp);


  sprintf(fname,"Q_mockmean.covar");
  fp = openfile(fname);
  for(j=1;j<=rsdq.np;++j)
    for(k=1;k<=rsdq.np;++k)
      fscanf(fp,"%lf",&rsdq.covar[j][k]);
  fclose(fp);

  // inver teh matrix
  tmp=dmatrix(1,rsdq.np,1,1);
  tmp2=dmatrix(1,rsdq.np,1,rsdq.np);
  for(i=1;i<=rsdq.np;++i)
    for(j=1;j<=rsdq.np;++j)
      tmp2[i][j]=rsdq.covar[i][j];
  gaussj(tmp2,rsdq.np,tmp,1);
  for(i=1;i<=rsdq.np;++i)
    for(j=1;j<=rsdq.np;++j)
      rsdq.covar[i][j]=tmp2[i][j];
  free_dmatrix(tmp,1,rsdq.np,1,1);
  free_dmatrix(tmp2,1,rsdq.np,1,rsdq.np);


}
