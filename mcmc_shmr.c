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
void mcmc_clustering_input(void);
double mcmc_shmr_initialize(double *a, double **cov1, double *avg1, double *start_dev);
double chi2_smf_shmr(void);
double chi2_wp_shmr(void);
int check_shmr_parameters(double *a);
int check_shmr_parameters_magnitude(double *a);
void input_shmr_chain(double *eval, double **evect);

/* external functions
 */
void covar_pca(int ii);
double external_constraints(void);
double chi2_smf_primus();

void mcmc_shmr()
{
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=400,NSTEP_MAX=1000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b;
 
  FILE *fpmc;
  char filename[1000];

  /* Get the SEED for the chain
   */
  IDUM=IDUM_MCMC;
  if(ARGC>2)
    IDUM = -1*abs(atoi(ARGV[2]));
  fprintf(stdout,"mcmc> ISEED= %d\n",(int)IDUM);
  
  /* CHeck to see if we need to reset the redshift
   */
  if(IDUM<=-1000)
    {
      i = abs(IDUM)%100;
      IZED = i;
      REDSHIFT = 0.2 + (i+0.5)*0.02;
      SIGMA_8 = SIGMA_8Z0*growthfactor(REDSHIFT);
      RESET_COSMOLOGY++;
      printf("WARNING: changing to z=%f and s8=%f\n",REDSHIFT,SIGMA_8);
      //if(i==0 || i==1 || i==2 || i==3 || i==5 || i==10) { wpl.stepfac_burn /= 10; }
    }
  

  /* what are the number of free parameters?
   * if doing no-color fits, params=11. If doing color fits, params=27.
   */
  n = 11 + 3*SATELLITE_PARAMETERIZATION + 2*VARIABLE_ALPHA;
  if(COLOR) n = 27 + 6*SATELLITE_PARAMETERIZATION + 4*VARIABLE_ALPHA;
  if(VARIABLE_EXCLUSION)n++;

  // if fitting BOSS cutoff, add two more params
  FITTING_BOSS = 1;
  if(FITTING_BOSS)n+=2;

  fprintf(stderr,"mcmc> [%d] free parameters.\n",n);
  wpl.ncf = n;
  MCMC_OUTPUT = 1;

  srand48(32498793);
  MCMC=Task.MCMC;

  pcheck=calloc(ptot,sizeof(int));
  if(DONT_FIT_CLUSTERING==0)mcmc_clustering_input();

  sprintf(filename,"%s_%d.mcmc",Task.root_filename,abs(IDUM));
  fpmc = fopen(filename,"w");


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

  chi2prev=mcmc_shmr_initialize(a,cov1,avg1,start_dev);
  niter++;

  for(i=1;i<=n;++i)
    {
      aprev[i] = a[i];
      chain[1][i] = a[i];
    }

  if(RESTART)
    {
      printf("mcmc> RESTART_FILE= [%s]\n",RESTART_FILE);
      input_shmr_chain(eval,evect);
      convergence = 1;
      goto SKIP_BURN;
    }

  pcnt=0;
  pcheck[pcnt]=1;

  stepfac=1;
  while(niter<NSTEP)
    {
      stepfac = wpl.stepfac_burn;
      for(i=1;i<=n;++i)
	wpl.a[i] = a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];
      wpl.reset_inversion = 1;

      if(check_shmr_parameters_magnitude(a))continue;

      if(MCMC_OUTPUT) {
	printf("START %d ",count+1);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("\n");
	fflush(stdout);
      }

      chi2a = chi2_wp_shmr();
      chi2b = chi2_smf_shmr();
      chi2 = chi2a + chi2b;

      if(MCMC_OUTPUT){
	printf("TRY %d ",++count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("%e\n",chi2);fflush(stdout);
      }

      pcheck[pcnt]=1;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  pcheck[pcnt]=0;
	  continue;	  
	}

      niter++;
      iweight[niter]++;

      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      //NB! BIAS FITTING ONLY
      chain[niter][6] = 0.25 + gasdev(&IDUM)*0.1;
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      fprintf(fpmc,"%d %d ",niter,count);
      for(i=1;i<=n;++i)
	fprintf(fpmc,"%e ",a[i]);
      fprintf(fpmc,"%e %e %e\n",chi2a,chi2b,chi2);fflush(fpmc);
      
    }

 SKIP_BURN:
  NSTEP = niter;
  stepfac=1.6/sqrt(n)*wpl.stepfac;
  while(niter<imax_chain)
    {

      if(convergence)goto SKIP_MATRIX;

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
      if(check_shmr_parameters_magnitude(a))continue;

      if(MCMC_OUTPUT) {
	printf("START %d ",count+1);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("\n");
	fflush(stdout);
      }


      chi2a = chi2_wp_shmr();
      chi2b = chi2_smf_shmr();
      chi2 = chi2a+chi2b;

      if(MCMC_OUTPUT){
	printf("TRY %d ",++count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("%e\n",chi2);fflush(stdout);
      }
      
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	continue;

      niter++;
      if(!convergence)NSTEP = niter;
      iweight[niter]++;

      if(niter==NSTEP_MAX)
	{
	  convergence = 1;
	  for(i=1;i<=n;++i)
	    eval_prev[i] = eval[i];
	}


      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      //NB! BIAS FITTING ONLY
      chain[niter][6] = 0.25 + gasdev(&IDUM)*0.1;
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      fprintf(fpmc,"%d %d ",niter,count);
      for(i=1;i<=n;++i)
	fprintf(fpmc,"%e ",a[i]);
      fprintf(fpmc,"%e %e %e\n",chi2a,chi2b,chi2);fflush(fpmc);
    }
}

double mcmc_shmr_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m;
  long IDUM = -556;
  FILE *fp;
  char aa[100];

  for(i=1;i<=wpl.ncf;++i)
    a[i] = wpl.a[i];

			 
  
  if(ARGC>3)
    {
      fp = openfile(ARGV[3]);
      fscanf(fp,"%d %d",&i,&j);
      for(i=1;i<=wpl.ncf;++i)
	fscanf(fp,"%lf",&a[i]);
      fclose(fp);
      for(i=1;i<=wpl.ncf;++i)wpl.a[i] = a[i];
    }
  

  start_dev[1]=0.03/10;                   // Mhalo norm
  start_dev[2]=0.03/10;                   // Mstellar norm
  start_dev[3]=0.03;                   // beta
  start_dev[4]=0.09;                   // delta
  start_dev[5]=0.3;                    // gamma
  start_dev[6]=0.03;     // sigma 
  start_dev[7] = 0.5;                // Mcut norm
  start_dev[8] = 2.0;                // M1 norm
  start_dev[9] = 0.3;                // Mcut slope
  start_dev[10] = 0.2;              // M1 slope
  start_dev[11] = 0.05; // 11p model: vary the satellite slope
  if(SATELLITE_PARAMETERIZATION)
    {
      start_dev[12] = 0.2;
      start_dev[13] = 0.1;
      start_dev[14] = 0.2;
    }
  if(VARIABLE_EXCLUSION) 
    {
      EXCLUSION_RADIUS = a[15];
      start_dev[15] = 0.1;
    }

  check_shmr_parameters_magnitude(wpl.a);

  if(!ThisTask)
    {
      printf("INITIAL VALUES: ");
      for(i=1;i<=wpl.ncf;++i)printf("%e ",a[i]);
      printf("\n");
    }

  for(i=1;i<=wp.ncf;++i)
    {
      avg1[i]=a[i];
      for(j=1;j<=wp.ncf;++j)
	cov1[i][j]=a[i]*a[j];
    }

  muh(0);
  x2=chi2_smf_shmr();
  x1=chi2_wp_shmr();

  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e\n",x1+x2);fflush(stdout);
    printf("INITIAL CHI2: %e %e\n",x1,x2);
    fflush(stdout);
  }
  //exit(0);
  return(x1+x2);
}

void mcmc_clustering_input()
{
  FILE *fp, *fp1, *fpc;
  int i,j,k,n;
  char filename[1000],aa[1000];
  double **covar, **tmp, **tmp2;

  if(DONT_FIT_CLUSTERING)return;

  // First, get the stellar mass bins
  fprintf(stdout,"mcmc> Opening file: [%s]\n",Files.StellarMassBinsClustering);
  fp = openfile(Files.StellarMassBinsClustering);
  wpl.nsamples = filesize(fp);
  
  fprintf(stderr,"mcmc> Reading %d stellar mass bins.\n",wpl.nsamples);
  for(i=1;i<=wpl.nsamples;++i)
    fscanf(fp,"%lf %lf",&wpl.mstar_wplo[0][i], &wpl.mstar_wphi[0][i]);
  fclose(fp);

  //Second, get the filenames and correlation functions
  fp = openfile(Files.ClusteringFilenames);
  if(filesize(fp) != wpl.nsamples) {
    fprintf(stderr,"mcmc> ERROR! wrong number of clustering files/bins!!\n");
    exit(0); }
  wpl.ndata = ivector(1,wpl.nsamples);
  for(i=1;i<=wpl.nsamples;++i)
    {
      fscanf(fp,"%s",filename);
      fp1 = openfile(filename);
      wpl.ndata[i] = filesize(fp1);
      wpl.rdata[i] = dvector(1,wpl.ndata[i]);
      wpl.xdata[i] = dvector(1,wpl.ndata[i]);
      wpl.edata[i] = dvector(1,wpl.ndata[i]);
      for(j=1;j<=wpl.ndata[i];++j)
	{
	  fscanf(fp1,"%lf %lf %lf",&wpl.rdata[i][j],&wpl.xdata[i][j],&wpl.edata[i][j]);
	  fgets(aa,1000,fp1);
	}
      fclose(fp1);
      fprintf(stderr,"mcmc> read [%d] lines from [%s]\n",wpl.ndata[i],filename);
    }

  //Thrid, get the covariance matrices if required
  if(!COVAR && !PCA)return;

  fp = openfile(Files.ClusteringCovarFilenames);
  if(filesize(fp) != wpl.nsamples) {
    fprintf(stderr,"mcmc> ERROR! wrong number of covar files/bins!!\n");
    exit(0); }
  for(i=1;i<=wpl.nsamples;++i)
    {
      fscanf(fp,"%s",filename);
      fp1 = openfile(filename);
      if(filesize(fp1)!=wpl.ndata[i]*wpl.ndata[i]) {
	fprintf(stderr,"mcmc> ERROR! wrong number elements in covar file [%s] %d %d\n",filename,filesize(fp1),wpl.ndata[i]);
	exit(0);
      }
      wpl.covar[i] = dmatrix(1,wpl.ndata[i],1,wpl.ndata[i]);
      for(j=1;j<=wpl.ndata[i];++j)
	for(k=1;k<=wpl.ndata[i];++k)
	  {
	    fscanf(fp1,"%lf",&wpl.covar[i][j][k]);
	    fgets(aa,1000,fp1);
	  }
      fclose(fp1);
      fprintf(stderr,"mcmc> read [%d] lines from [%s]\n",wpl.ndata[i]*wpl.ndata[i],filename);
      covar_pca(i);
    }

  // Last, invert the covariance matrices
  fprintf(stderr,"mcmc> inverting covariance matrices.\n");

  for(i=1;i<=wpl.nsamples;++i)
    {
      n = wpl.ndata[i];
      tmp=dmatrix(1,n,1,1);
      tmp2=dmatrix(1,n,1,n);
      for(j=1;j<=n;++j)
	for(k=1;k<=n;++k)
	  tmp2[j][k]=wpl.covar[i][j][k];
      gaussj(tmp2,n,tmp,1);
      for(j=1;j<=n;++j)
	for(k=1;k<=n;++k)
	  wpl.covar[i][j][k]=tmp2[j][k];
      free_dmatrix(tmp,1,n,1,1);
      free_dmatrix(tmp2,1,n,1,n);
    }

}


double chi2_wp_shmr()
{
  static niter=0;
  double xm[1000];
  int i,j,k;
  double chi2,chi2tot=0,mlo,mhi;
  FILE *fp;

  if(DONT_FIT_CLUSTERING)return 0;

  for(i=1;i<=wpl.nsamples;++i)
    {
      mlo = pow(10.0,wpl.mstar_wplo[0][i]);
      mhi = pow(10.0,wpl.mstar_wphi[0][i]);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);

      GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      BETA = pow(OMEGA_Z,GAMMA)/GALAXY_BIAS;

      RESET_FLAG_1H = RESET_FLAG_2H = 1;
      RESET_KAISER++;

      chi2 = 0;
      for(j=1;j<=wpl.ndata[i];++j)
	{
	  xm[j] = projected_xi(wpl.rdata[i][j]);
	  chi2 += (xm[j]-wpl.xdata[i][j])*(xm[j]-wpl.xdata[i][j])/(wpl.edata[i][j]*wpl.edata[i][j]);
	  if(MCMC_OUTPUT)
	    printf("WP%d %d %e %e %e %e\n",niter,i,wpl.rdata[i][j],wpl.xdata[i][j],wpl.edata[i][j],xm[j]);
	}
      if(COVAR && !PCA)
	{
	  chi2 = 0;
	  for(j=1;j<=wpl.ndata[i];++j)
	    for(k=1;k<=wpl.ndata[i];++k) {
	      chi2 += (wpl.xdata[i][j]-xm[j])*(wpl.xdata[i][k]-xm[k])*wpl.covar[i][j][k];
	      //printf("BOO %d %d %e %e %e\n",j,k,wpl.covar[i][j][k],(wpl.xdata[i][j]-xm[j])*(wpl.xdata[i][k]-xm[k])*wpl.covar[i][j][k], chi2);
	    }
	}
      if(PCA)
	chi2 = pca_chi2(xm,i);
      if(MCMC_OUTPUT)
	printf("WPCHI%d %d %e\n",niter,i,chi2);
      chi2tot += chi2;
    }
  if(EXTERNAL_CONSTRAINTS)
    chi2tot += external_constraints();

  niter++;
  return chi2tot;
}


double chi2_smf_shmr()
{
  static int n=0, npr, npb, n_ratio, nr, k, istart_b, istart_r;
  static double *mstar, *nstar, *estar, **covar, *eval, **evect;
  static double *mstar_b, *nstar_b, *estar_b, **covar_b;
  static double *mstar_r, *nstar_r, *estar_r, **covar_r;
  static double *mstar_ratio, *nstar_ratio, *estar_ratio, **covar_ratio;
  int i, j, i1, ii, ibuf, nlim = 20, eflag;
  FILE *fp;
  static double dlog10m;
  static int niter=0;
  double ng1, ng2, mlo, mhi, ngal, ngalprev, dmsdmh, m1, m2, ns1, ns2, **tmp, **tmp2, mass, nsat, chi1, cbias;
  char fname[1000], aa[1000];
  double chi2ngal, xmodel[1000], xmodelr[100],xmodelb[100];
  FILE *outfile, *outfile_cov;

  if(DONT_FIT_SMF==2)return chi2_smf_primus();

  /* The HOD set up is different for the stellar mass function because
   * only really need to know the total number of central and satellites
   * not the full Ncen=f(M) and Nsat=f(sat) (doesn't call these functions)
   */
  if(!n)
    {
      fp = openfile(Files.StellarMassFunctionFilename);
      n = filesize(fp);
      nlim = n;
      mstar = dvector(1,nlim);
      nstar = dvector(1,nlim);
      estar = dvector(1,nlim);
      
      for(i=1;i<=n;++i) {
	fscanf(fp,"%lf %lf %lf",&mstar[i],&nstar[i],&estar[i]);
	fgets(aa,1000,fp); }
      fclose(fp);
      fprintf(stderr,"Read [%d] lines from [%s]\n",n,Files.StellarMassFunctionFilename);
	  
      if(COVAR>1)
	{	
	  fp = openfile(Files.StellarMassFunctionCovarFilename);
	  if(filesize(fp)!=n*n)
	    {
	      fprintf(stderr,"FILESIZE mismatch: cij and data not match in SMF: %d %d\n",n*n,filesize(fp));
	      exit(0);
	    }
	  covar = dmatrix(1,nlim,1,nlim);
	  for(i=1;i<=n;++i)
	    for(j=1;j<=n;++j)
	      fscanf(fp,"%lf",&covar[i][j]);
	  fclose(fp);	
	  fprintf(stderr,"Read [%d] lines from [%s]\n",n*n,Files.StellarMassFunctionCovarFilename);
	  
	  printf("INVERTING COVARIANCE MATRIX %d\n",n);
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
	}
     
      // using the tabulated data, determine the binsize
      dlog10m = mstar[2] - mstar[1];
      fprintf(stderr,"mcmc> SMF binsize: %f\n",dlog10m);

      if(SMF_PCA)
	{
	  fp = openfile(Files.StellarMassFunctionCovarFilename);
	  if(filesize(fp)!=n*n)
	    {
	      fprintf(stderr,"FILESIZE mismatch: cij and data not match in SMF: %d %d\n",n*n,filesize(fp));
	      exit(0);
	    }
	  covar = dmatrix(1,nlim,1,nlim);
	  eval = dvector(1,nlim);
	  evect = dmatrix(1,nlim,1,nlim);
	  for(i=1;i<=n;++i)
	    for(j=1;j<=n;++j)
	      fscanf(fp,"%lf",&covar[i][j]);
	  fclose(fp);	
	  fprintf(stderr,"Read [%d] lines from [%s]\n",n*n,Files.StellarMassFunctionCovarFilename);

	  for(i=1;i<=n;++i)
	    for(j=1;j<=n;++j)
	      covar[i][j] /= sqrt(covar[i][i]*covar[j][j]);

	  /* Do the singular-value decomposition on the covariance matrix
	   */
	  svdcmp(covar,n,n,eval,evect);
	}

    }

  chi2ngal = 0;
  ngal = 0;

  wpl.reset_inversion = 1;
  i=set_up_hod_for_shmr(1.0E+10+drand48()*1.0E+10,3.0E10,wpl.a);
  //if we have an issue, return big chi2 value
  if(i==-1){ niter++; return 1.0E+7; }
  
  for(i=1;i<=n;++i)
    xmodel[i] = 0;

  for(i=1;i<=n;++i)
    {
      mlo = pow(10.0,mstar[i]-dlog10m/2); // these are the bin limits (don't use wpl here since Ncen is not called)
      mhi = pow(10.0,mstar[i]+dlog10m/2);
      
      eflag = ngal = cbias = 0;
      if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
      //ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      
      //NB!! BIAS FITTING ONLY (commented out above ngal)
      if(!eflag)
      	{
	  // NB! I've also changed the halo mass function to account for satellites.
	  ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;      
	}

      ngalprev = ngal;
      ngal = ngal*pow(HUBBLE_UNITS,3.0)/dlog10m; // Galaxy density, convert to h72 units
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
	  
      if(COVAR<2)
	chi2ngal += (ngal-nstar[i])*(ngal-nstar[i])/(estar[i]*estar[i]);
      xmodel[i] = ngal;
      
      if(MCMC_OUTPUT)
	printf("SMF%d %e %e %e %e %e\n",niter,mstar[i],nstar[i],estar[i],ngal,ms_to_mhalo(pow(10.0,mstar[i]),wpl.a));

      if(ERROR_FLAG == 1){
	printf("\n STOPPED IN chi2_stellar_mass_function %d\n",i);
	printf("SMF%d %e %e %e %e %e\n",niter,mstar[i],nstar[i],estar[i],ngal,(ns1-ns2)*pow(HUBBLE_UNITS,3.0)/dlog10m);
	printf("%e %e %e %e\n",HOD.M_low, HOD.M_min, HOD.M_max, N_sat(HOD.M_max));
	ERROR_FLAG=0;
	//exit(0);
      }
      if(ngal<1.0E-15)break;
      
    }
      
  // now calculate the chi^2 with the covariance matrix
  if(COVAR==2)
    for(i=1;i<=n;++i)
      for(j=1;j<=n;++j){
	chi2ngal += (xmodel[i]-nstar[i])*(xmodel[j]-nstar[j])*covar[i][j];
      }

  if(SMF_PCA)
    {
      chi2ngal = 0;
      for(i=1;i<=SMF_PCA;++i)
	{
	  chi1 = 0;
	  for(j=1;j<=n;++j)
	    chi1 += evect[j][i]*(xmodel[j]-nstar[j])/estar[j];
	  chi2ngal += chi1*chi1/eval[i];
	}
    }


  muh(0);
  // Set this flag to 1 to avoid fitting the SMF
  // This must go and end and not beginning because otherwise galaxy_density is not calcualted properly.
  // ### why was code failing when galaxy density not calculated ??? 
  if(DONT_FIT_SMF)
    chi2ngal=0;

  if(MCMC_OUTPUT)
    printf("SMFCHI%d %e\n",niter,chi2ngal);
  niter++;

  return chi2ngal;
}

/* check and make sure the parameters selected are within the 
 * acceptable parameter space.
 */
int check_shmr_parameters(double *a)
{
  int check_flag, ibuf;
  float m1;

  // do some things that need to be done each time
  if(VARIABLE_EXCLUSION) {
    EXCLUSION_RADIUS = a[wpl.ncf]; // the exclusion radius is always the last parameter
    if(EXCLUSION_RADIUS>2.0)return wpl.ncf;
    if(EXCLUSION_RADIUS<1.0)return -wpl.ncf;
  }
  wpl.reset_fred = 1;

  if(a[1]<=12.0)return 1;      // Mhalo_norm
  if(a[1]>=13.8)return -1;
    
  if(a[2]<=10.3)return 2;      // Mstellar norm
  if(a[2]>=12.0)return -2;
 
  if(a[3]<=0.2)return 3;       // Beta 
  if(a[3]>=0.8)return -3;

  if(a[4]<=0.01)return 4;       // delta (If <0, affects the threshold HOD, strange things happen .. don't go below 0.3)
  if(a[4]>=1.4)return -4;
  
  if(a[5]<=-0.5)return 5;      // gamma (can be <0)
  if(a[5]>=5.0)return -5;

  if(a[6]<=0.01)return 6;  // sigma_log_sm
  if(a[6]>=0.5)return -6;
  
  // PUTTING A PRIOR ON THIS BASED ON THE MEASUREMENT ERROR
  //if(a[6]<0.16)return 6;

  // FOR BIAS FITTING ONLY!!! NB NB NB (commenting out above if(a[6]...)
  //wpl.a[6] = a[6] = HOD.sigma_logM; //(putting it in a global variable for convenience)

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

  if(a[11]>3)return 11;
  if(a[11]<0.1)return -11;

  if(SHMR_PARAMS>11)
    {
      if(a[12]>1.0E17)return 12;
      if(a[12]<1.0E10)return -12;
      if(a[13]>3)return 13;
      if(a[13]<0)return -13;
      if(a[14]>1.0E17)return 14;
      if(a[14]<1.0E3)return -14;
    }

  return 0;
}

/*-----------------------------------------------------------------------------
 * 
 * INPUT AN EXTERNAL CHAIN FOR COVARIANCE MATRIX
 *
 *----------------------------------------------------------------------------*/

void input_shmr_chain(double *eval, double **evect)
{
  int np,i,j,n,total_weight,nrot,i1,k;
  double **chain, **cov1, *avg1, **tmp, **tmp1;
  FILE *fp;
  char aa[1000];

  printf("%s\n",RESTART_FILE);
  fp = openfile(RESTART_FILE);
  np = filesize(fp);
  chain = dmatrix(1,np,1,wpl.ncf);

  n = wpl.ncf;
  cov1=dmatrix(1,n,1,n);
  avg1=dvector(1,n);

  tmp=dmatrix(1,n,1,n);
  tmp1=dmatrix(1,n,1,1);

  for(i=1;i<=np;++i)
    {
      fscanf(fp,"%d %d",&i1,&i1);
      for(j=1;j<=wpl.ncf;++j)
	fscanf(fp,"%lf",&chain[i][j]);
      fgets(aa,1000,fp);
    }
  fclose(fp);

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


/* check and make sure the parameters selected are within the 
 * acceptable parameter space.
 */
int check_shmr_parameters_magnitude(double *a)
{
  int check_flag, ibuf;
  float m1;

  // do some things that need to be done each time
  if(VARIABLE_EXCLUSION) {
    EXCLUSION_RADIUS = a[wpl.ncf]; // the exclusion radius is always the last parameter
    if(EXCLUSION_RADIUS>2.0)return wpl.ncf;
    if(EXCLUSION_RADIUS<1.0)return -wpl.ncf;
  }
  wpl.reset_fred = 1;

  if(a[1]<=12.0)return 1;      // Mhalo_norm
  if(a[1]>=15.8)return -1;
    
  if(a[2]<=10.3)return 2;      // Mstellar norm
  if(a[2]>=13.0)return -2;
 
  if(a[3]<=0.2)return 3;       // Beta 
  if(a[3]>=0.8)return -3;

  if(a[4]<=0.01)return 4;       // delta (If <0, affects the threshold HOD, strange things happen .. don't go below 0.3)
  if(a[4]>=1.4)return -4;
  
  if(a[5]<=-0.5)return 5;      // gamma (can be <0)
  if(a[5]>=5.0)return -5;

  if(a[6]<=0.01)return 6;  // sigma_log_sm
  if(a[6]>=0.5)return -6;
  
  // PUTTING A PRIOR ON THIS BASED ON THE MEASUREMENT ERROR
  //if(a[6]<0.16)return 6;

  // FOR BIAS FITTING ONLY!!! NB NB NB (commenting out above if(a[6]...)
  //wpl.a[6] = a[6] = HOD.sigma_logM; //(putting it in a global variable for convenience)

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

  if(a[11]>3)return 11;
  if(a[11]<0.1)return -11;

  if(SHMR_PARAMS>11)
    {
      if(a[12]>1.0E17)return 12;
      if(a[12]<1.0E10)return -12;
      if(a[14]>3)return 13;
      if(a[14]<0)return -13;
      if(a[13]>1.0E17)return 14;
      if(a[13]<1.0E3)return -14;
    }

  return 0;
}
