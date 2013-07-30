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

/* Internal functions.
 */
double chi2_wp_wrapper(double *a);
void mcmc_restart2(double *start_dev, int np);
int mcmc_restart3(double **chain, int n, double *chi2_prev, int *iweight);

int USE_IWEIGHT = 0;

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
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=10000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b;

  opar=dvector(1,100);

  MCMC=Task.MCMC;

  pcheck=calloc(ptot,sizeof(int));

  if(MCMC>1 && !COVAR)
    wp.esys=0.05;

  if(!ThisTask)printf("ESYS %f %d\n",wp.esys,MCMC);
  wp_input();

  Work.imodel=2;
  Work.chi2=1;

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

  IDUM=IDUM_MCMC;

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

  stepfac=1;
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
      /* stepfac=0.7; */
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      
      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  /* if(HOD.free[++i])SIGV    = a[++j]; */
	}
      if(VBIAS_C<0)continue;

      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[3];

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
	printf("%e\n",chi2);fflush(stdout);
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

      if(!ThisTask){
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e\n",chi2);fflush(stdout);
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
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  /* if(HOD.free[++i])SIGV    = a[++j]; */
	}
      if(VBIAS_C<0)continue;
      
      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[3];

      /* Draw random value of cvir from prior.
       */
      /* CVIR_FAC = a[n]; */
      /* if(CVIR_FAC<0.3 || CVIR_FAC>1.2)continue; */
      /* CVIR_FAC = 0.7*drand48()+0.3; */
      /* GAMMA = gasdev(&IDUM)*0.02 + 0.15; */
      //      printf("GAMMA %d %f %f\n",count+1,GAMMA,CVIR_FAC);

      chi2a=chi2_wp_wrapper(a);
      chi2 = chi2a+chi2b;

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

      if(!ThisTask) {
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e\n",chi2);fflush(stdout);
	
	if(MCMC==1)
	  {
	    printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
		   number_weighted_central_mass(),
		   qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	    
	    fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	    printf("FSAT %d %e %e %e %e\n",niter,fsat,HOD.M_min,HOD.sigma_logM,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	  }
      }

    }
}

double chi2_wp_wrapper(double *a)
{
  static int flag=1;
  static double *b;
  int i,j;

  if(flag)
    {
      b=dvector(1,100);
      flag=0;
    }

  for(j=0,i=1;i<=N_HOD_PARAMS;++i) {
    if(HOD.free[i] && i!=5) { 
      if(a[++j]<=0) { printf("NEG %d %d %e\n",i,j,a[j]); return(1.0E7); } }
    if(HOD.free[i] && i==5) {
      ++j; }
  }

  i=0;j=0;
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_min */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M1 */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_cut */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* sigma_logM */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* cvir_fac */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* MaxCen */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_sat_break */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha1 */

  return(chi2_wp(b));
}

double mcmc_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m;
  long IDUM = -556;

  omega_m = 1;
  if(MCMC>1)
    omega_m = OMEGA_M;

  i=0;j=0;
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_min/omega_m);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M1/omega_m);start_dev[j]=0.001; } //.0005
  if(HOD.free[++i]){ a[++j]=HOD.alpha;start_dev[j]=0.03; } //.005
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_cut/omega_m);start_dev[j]=0.01; } //.001
  if(HOD.free[++i]){ a[++j]=log10(HOD.sigma_logM);start_dev[j]=0.01; }
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
      if(HOD.free[++i])a[++j]=SIGMA_8;
      if(HOD.free[++i])a[++j]=VBIAS;
      if(HOD.free[++i])a[++j]=VBIAS_C;
      if(HOD.free[++i])a[++j]=GAMMA;
      if(HOD.free[++i])a[++j]=SPECTRAL_INDX;
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
      if(HOD.free[++i]){ SIGMA_8 = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS   = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS_C = a[++j]; start_dev[j] = 0.02; } 
      if(HOD.free[++i]){ GAMMA   = a[++j]; start_dev[j] = 0.015; } 
      if(HOD.free[++i]){ SPECTRAL_INDX    = a[++j]; start_dev[j] = 0.02; }
    }

  x1=chi2_wp_wrapper(a);
  
  x2=0;

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
