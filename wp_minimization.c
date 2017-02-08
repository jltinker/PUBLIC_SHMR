#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif


#include "header.h"

/* This is a series of routines to calculate the projected correlation
 * function from the real-space one (HOD calculation) and then chi2
 * minimize the parameters based on the SDSS data.
 */
int FIT_WITH_COSMO = 0,
  FIRST_CALL = 0;

double rp_g1;
double func_wp(double z);
double func_wp_1h(double z);
double func_wp_rspace(double z);
double func_wp_matter(double z);
double chi2_wp(double *a);

void wp_input(void);
void initial_wp_values(double *a, double **pp, double *yy);

/* The two functions below calculate the projected correlation
 * function from the HOD real-space correlation function at 
 * the input value of r_p.
 */

/* This integrates the calculated real-space correlation function
 * along the line-of-sight to get the projected correlation function (wp(rp).
 * The value passed if the projected separation.
 *
 * In most SDSS work, the maximum line of sight separation considered when
 * claculating wp is pi=40 Mpc/h, which is the default, but can be set
 * in the parameter file if desired.
 */
double projected_xi(double r)
{
  double x,zmax;

  rp_g1=r*r;
  two_halo_real_space(1.0);
  zmax=sqrt(XI_MAX_RADIUS*XI_MAX_RADIUS-rp_g1);
  if(zmax>wp.pi_max)zmax=wp.pi_max;
  zmax = wp.pi_max;
  x=2*qromo(func_wp,log(0.001),log(zmax),midpnt);
  /* x=2*qtrap(func_wp,log(0.001),log(zmax),1.0E-4); */
  return(x);
}

double projected_xi_1halo(double r)
{
  double x,zmax;

  rp_g1=r*r;
  two_halo_real_space(1.0);
  zmax=sqrt(XI_MAX_RADIUS*XI_MAX_RADIUS-rp_g1);
  if(zmax>wp.pi_max)zmax=wp.pi_max;
  zmax = wp.pi_max;
  x=2*qromo(func_wp_1h,log(0.001),log(zmax),midpnt);
  return(x);
}

double projected_xi_rspace(double r)
{
  double x,zmax;

  rp_g1=r*r;
  two_halo_real_space(1.0);
  zmax=sqrt(XI_MAX_RADIUS*XI_MAX_RADIUS-rp_g1);
  if(zmax>wp.pi_max)zmax=wp.pi_max;
  zmax = wp.pi_max;
  //x=2*qromo(func_wp_rspace,log(0.001),log(zmax),midpnt);
   x=2*qtrap(func_wp_rspace,log(0.001),log(zmax),1.0E-4); 
  return(x);
}

/* Same as above but for matter.
 * I've set the maximum line-of-sight separation to be 50 Mpc/h,
 * which should be good for visual comparison purposes.
 */
double projected_xi_matter(double r)
{
  double x,zmax;

  rp_g1=r*r;
  zmax=50.0;
  x=2*qtrap(func_wp_matter,log(0.001),log(zmax),1.0E-3);
  return(x);
}

/* Function called from qromo/qtrap to get xi->wp
 */
double func_wp(double z)
{
  double r;
  z=exp(z);
  r=sqrt(rp_g1 + z*z);

  return(z*(one_halo_real_space(r)+linear_kaiser_distortion(r,z)));
}

/* Function called from qromo/qtrap to get xi->wp
 */
double func_wp_1h(double z)
{
  double r;
  z=exp(z);
  r=sqrt(rp_g1 + z*z);
  return(z*(one_halo_real_space(r)));
}

/* Function called from qromo/qtrap to get xi->wpm but without 
 * correction for redshift-space distortions in the two-halo term
 */
double func_wp_rspace(double z)
{
  double r;
  z=exp(z);
  r=sqrt(rp_g1 + z*z);
  return(z*(one_halo_real_space(r)+two_halo_real_space(r)));
}
  
/* Function called from qromo/qtrap to get xi->wp (dark matter)
 */
double func_wp_matter(double z)
{
  double r;
  z=exp(z);
  r=sqrt(rp_g1 + z*z);
  return(z*(xi_interp(r)));
  return(z*linear_kaiser_distortion(r,z));
}
  
/********************************************************************
 * Below is the actual minimization of the wp data.
 * Relevant variables: [default]
 * 
 * COVAR ->  [1]=use covariance matrixl; 0=diagonal error bars
 * DEPROJECTED -> 1=we're fitting a real-space xi(r) [0]= fitting wp(rp)
 *
 * HOD.free[] is a vector which holds 1/0 as to whether or not a parameter is going
 * to be help constant during the chi^2 minimization. 1==vary 0==constant. The default
 * on any will be [0].
 *
 *  i     variable
 * ---    --------
 * [1] ->  M_min
 * [2] ->  M1
 * [3] ->  alpha
 * [4] ->  M_cut
 * [5] ->  sigmaM
 * [6] ->  CVIR_FAC
 * [7] ->  MaxCen (M_cen_max)
 *
 * Currently I have no checks on whether the values of this vector line
 * up correctly with the specific HOD pdfs, so double-check hod.bat files.
 *
 * Once the code is finished, it will output the values of the HOD parameters
 * to a file called [filename].fit (+ the bias, satellite fraction, and chi^2).
 * Then it outputs the mean <N>_M to a file called [filename].hod.
 */

void wp_minimization(char *fname)
{
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m;
  FILE *fp;
  char aa[1000];

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF W_P(R_P) DATA..........\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  //OUTPUT = 0;
  FIRST_CALL = 1;

  if(POWELL)
    FTOL=1.0E-3;
  else
    FTOL=1.0E-4;

  for(n=0,i=1;i<=N_HOD_PARAMS;++i)
    {
      n+=HOD.free[i];
      if(!OUTPUT)continue;
      printf("wp_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  if(XCORR)n*=2;
  if(OUTPUT)printf("wp_min> Number of free parameters: %d\n",n);

  printf("FNAME %s\n",Task.root_filename);
  wp_input();

  wp.ncf=n;
  a=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);


  initial_wp_values(a,pp,yy);
  printf("IVALS %e %e %e %e\n",a[1],a[2],a[3],a[4]);


  if(POWELL) 
    {
      if(OUTPUT)printf("wp_min> starting powell.\n");
      powell(a,pp,n,FTOL,&niter,&chi2min,chi2_wp);
      chi2min = chi2_wp(a);
    }
  else
    {
      if(OUTPUT)printf("wp_min> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_wp,&niter);
      for(i=1;i<=n;++i)a[i]=pp[1][i];
      chi2min = chi2_wp(a);
    }	

  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;

  printf("POWELL %e %e ",chi2min,HOD.M_min);
  for(i=1;i<=n;++i)printf("%e ",a[i]);
  printf(" %f\n",GALAXY_BIAS);

  /* These outputs are for easy cut & paste into
   * another batch file.
   */
  //output_parameter_file(fname);

  /* Output the fit and the HOD curve.
   */
  printf("FNAME2 %s\n",Task.root_filename);
  sprintf(aa,"%s.fit",Task.root_filename);
  fp=fopen(aa,"w");
  fprintf(fp,"%e %e ",chi2min,HOD.M_min);
  for(i=1;i<=n;++i)fprintf(fp,"%e ",a[i]);
  fprintf(fp," %f\n",GALAXY_BIAS);
  fclose(fp);

  sprintf(aa,"%s.HOD",Task.root_filename);
  fp=fopen(aa,"w");
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      fprintf(fp,"%e %e %e %e\n",m,N_cen(m),N_sat(m),N_avg(m));
    }
  fclose(fp);

  fprintf(stderr,"here\n");
  free_dvector(a,1,n);
  if(POWELL)
    free_dmatrix(pp,1,n,1,n);
  else
    free_dmatrix(pp,1,n+1,1,n);
  free_dvector(yy,1,n+1);
  fprintf(stderr,"here\n");

  free_dvector(wp.r,1,wp.np);
  free_dvector(wp.x,1,wp.np);
  free_dvector(wp.e,1,wp.np);
  if(COVAR)
    free_dmatrix(wp.covar,1,wp.np,1,wp.np);
  fprintf(stderr,"done in wp_min\n");
}

double integrated_wp_bin(double r)
{
  return(projected_xi(r)*r);
}
double integrated_xi_bin(double r)
{
  //return((1+one_halo_real_space(r))*r*r);
  return((one_halo_real_space(r)+two_halo_real_space(r))*r*r);
}

double chi2_wp(double *a)
{
  static int flag=1,niter=0,ichi=-1;
  static double x[100],mmin_prev=0,t0=-1,t1,sig_prev=0,chi2_prev,chi2_array[10];
  double **tmp,**tmp2,chi2,x1,ta1,ta2,dt1h,dt2h,par_chi,chi2ngal;
  int i,j,k,ncf_hod;

  double rlo,rhi,rmin,rmax,dlogr,integrated_wp_bin(),omega_temp;


  if(FIRST_CALL)
    {
      flag = 1;
      FIRST_CALL = 0;
    }

  t0 = clock();
  if(HOD.free[1])FIX_PARAM = 0;

  wp.iter=niter;

  for(j=0,i=1;i<=N_HOD_PARAMS;++i)
    if(HOD.free[i])
      if(a[++j]<=0) { printf("%d %e\n",j,a[j]); return(1.0E7); }
  ncf_hod = j;

  RESET_FLAG_1H=1;
  RESET_FLAG_2H=1;
  RESET_KAISER++;

  i=0;j=0;
  if(HOD.free[++i])HOD.M_min=a[++j];
  if(HOD.free[++i])HOD.M1=a[++j];
  if(HOD.free[++i])HOD.alpha=a[++j];
  if(HOD.free[++i])HOD.M_cut=a[++j];
  if(HOD.free[++i])HOD.sigma_logM=a[++j];


  /* NORMALIZE the halo masses by OMEGA_M
   */
  /*
  if(MCMC>1)
    {
      HOD.M1 *= OMEGA_M;
      HOD.M_cut *= OMEGA_M;
    }
  */
  /* if(HOD.i_wp<6 && !POWELL)HOD.sigma_logM = 0.15; */

  /* HUBBLE BUBBLE-- only for low-L samples
   */
  //if(HOD.sigma_logM>0.5 && GALAXY_DENSITY>2.0e-3)return(1.1e7);

  /* UBER HOD: for -19.5 and -20.0 have sigma_logM<0.25)
   *
  if(HOD.pdfc!=1 && HOD.pdfc != 7) {
    if(HOD.sigma_logM>0.25 && GALAXY_DENSITY>0.0030)return(1.1e7);
    if(HOD.sigma_logM<0.15 && GALAXY_DENSITY>0.0030)return(1.1e7);
  }
  */

  if(HOD.pdfc!=9) {
    if(HOD.free[++i])CVIR_FAC=a[++j];
    if(HOD.pdfc>=7) {
      if(HOD.free[++i])HOD.M_cen_max=a[++j]; }
    else {
      if(HOD.free[++i])HOD.MaxCen=a[++j]; }
  }
  if(HOD.free[++i])HOD.M_sat_break=a[++j];
  if(HOD.free[++i])HOD.alpha1=a[++j];
  

  if(!ThisTask) {
    printf("START %d ",niter);
    for(i=1;i<=ncf_hod;++i)printf("%e ",a[i]);
    printf("\n");
  }
  

  if(HOD.pdfs==2 && HOD.M_cut<1.0e7)return(1.0e7);
  if(HOD.pdfs==2 && HOD.M_cut>1.0e15)return(1.0e7);

  /* if(HOD.M_min>HOD.M_max)return(1.0e7); */

  /* I've noticed some problems when sigma_logM gets to be
   * unresonably high or low, so I've put some limits on the
   * values they can have when doing mag-bin fitting.
   */
  if(HOD.pdfc==6) {
    if(HOD.sigma_logM<0.07)return(1.0e7);
    if(HOD.sigma_logM>1.2)return(1.0e7);
  }
  if(HOD.pdfc==2 || HOD.pdfc==9) {
    if(HOD.sigma_logM>1.8)return(1.0e7);
    if(HOD.sigma_logM<0.05)return(1.0e7);
  }
  if(HOD.M1>1.0e17)return(1.0e7);

  if(FIX_PARAM==2)
    {
      HOD.M1=HOD.M_low;
      x1=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      if(x1<GALAXY_DENSITY)return(1.0e7);
      HOD.M1=pow(10.0,14.8);
      x1=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      if(x1>GALAXY_DENSITY)return(1.0e7);
      HOD.M1=0;
    }


  /* Check the make sure these are reasonable parameters
   * (Assuming that M_min is NOT a FREE parameter but is
   * calculated from the GALAXY_DENSITY.)
   */
  if(FIX_PARAM==1 && !HOD.color)
    {
      HOD.M_min=pow(10.0,8.0);
      HOD.M_low=set_low_mass();
      if(HOD.M_low<1.0e8)HOD.M_low=1.0e8;
      x1=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      //      fprintf(stderr,"PC1 %e %e\n",x1,GALAXY_DENSITY); 
      if(x1<GALAXY_DENSITY) {
	fprintf(stdout,"PCHECK %e %e %e\n",x1,GALAXY_DENSITY,HOD.M_low);
	return(1.0e7); }
      HOD.M_min=pow(10.0,14.8);
      if(HOD.pdfc==7 && HOD.pdfc==8)
	HOD.M_min=HOD.M_cen_max*0.99; 
      HOD.M_low=set_low_mass();
      x1=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      /* fprintf(stderr,"PC2 %e %e %e %e\n",HOD.M_min,HOD.M_low,x1,GALAXY_DENSITY); */
      if(x1>GALAXY_DENSITY) { 
	fprintf(stdout,"PCHECK %e %e\n",x1,GALAXY_DENSITY);
	return(1.0e7); }
      HOD.M_min=0;
    }

  if(ERROR_FLAG)
    {
      ERROR_FLAG=0;
      return(1e7);
    }

  if(HOD.free[0] || HOD.free[1])
    GALAXY_DENSITY=0;

  set_HOD_params();

  if(HOD.free[0])
    {
      chi2ngal = (GALAXY_DENSITY-wp.ngal)*(GALAXY_DENSITY-wp.ngal)/wp.ngal_err/wp.ngal_err;
      if(chi2ngal>1.0E3)return(chi2ngal);
    }

  if(ERROR_FLAG)
    {
      ERROR_FLAG=0;
      return(1e7);
    }

  /* if(HOD.pdfs==3 && HOD.M_cut<HOD.M_low)return(1.1e7); */
  if(HOD.M_min>HOD.M1)return(1.0e7);
  /* if(HOD.pdfs==3)if(HOD.M_cut<HOD.M_min*0.9)return(1.2e7); */

  mmin_prev=HOD.M_min;
  sig_prev=HOD.sigma_logM;

  
  if(HOD.color>=niter)
    flag = 1;
  if(MCMC>2)
    flag = 0; //we've already onverted the matrix somewhere else

  if(flag && COVAR)
    {
      printf("INVERTING COVARIANCE MATRIX\n");
      flag=0;
      tmp=dmatrix(1,wp.np,1,1);
      tmp2=dmatrix(1,wp.np,1,wp.np);
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  tmp2[i][j]=wp.covar[i][j];
      gaussj(tmp2,wp.np,tmp,1);
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  wp.covar[i][j]=tmp2[i][j];
      free_dmatrix(tmp,1,wp.np,1,1);
      free_dmatrix(tmp2,1,wp.np,1,wp.np);
      
    }
 
  rmax = wp.r[wp.np];
  rmin = wp.r[1];
  dlogr = (log(rmax)-log(rmin))/(wp.np-1);

  printf("%e %e %e\n",one_halo_real_space(1.0),two_halo_real_space(1.0),projected_xi(1.0));
  x[1] = projected_xi(1.0);
  printf("%e\n",x[1]);

  omega_temp = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  //omega_temp = OMEGA_M;
  BETA = pow(omega_temp,GAMMA)/qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)*
    GALAXY_DENSITY;
  if(OUTPUT)
    printf("BETA = %f %f\n",BETA,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/
    GALAXY_DENSITY);

  rlo = exp(log(rmin) - 0.5*dlogr);
  for(i=1;i<=wp.np;++i)
    {
      rhi = exp(dlogr)*rlo;
      switch(DEPROJECTED) {
      case 1:
	x[i]=one_halo_real_space(wp.r[i])+two_halo_real_space(wp.r[i]);
	break;
      case 2:
	x[i] = wtheta(wp.r[i]);
	break;
      case 3: case 0:
	x[i]=projected_xi(wp.r[i]);      
	if(wp.format==3)
	  x[i]/=wp.r[i];
	break;
      }
      if(OUTPUT && !ThisTask)
	printf("WP%d %f %e %e %e %e\n",niter,wp.r[i],wp.x[i],x[i],rlo,rhi);

      rlo=rhi;

    }
  if(ERROR_FLAG && DENSITY_DEPENDENCE)
    ERROR_FLAG=0;

  if(ERROR_FLAG)
    {
      ERROR_FLAG=0;
      return(1e7);
    }

  chi2=0;

  if(COVAR)
    {
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  {
	    x1=(x[i]-wp.x[i])*(x[j]-wp.x[j])*wp.covar[i][j];
	    chi2+=x1;
	    /*
	    if(OUTPUT==2 && !ThisTask)
	      printf("CHI%d %e %e %e %e %e %e %e\n",
		     niter+1,wp.x[i],x[i],wp.x[j],x[j],wp.covar[i][j],x1,chi2);
	    if(MCMC>1 && i==j && !ThisTask)
	      printf("CHIWP %d %e %e %e %e\n",
		     niter+1,wp.r[i],x[i],wp.x[i],wp.e[j]);
	    */
	  }
    }
  else
    {
      if(!PCA) 
	{
	  for(i=1;i<=wp.np;++i)
	    {
	      x1=(x[i]-wp.x[i])*(x[i]-wp.x[i])/
		(wp.e[i]*wp.e[i] + wp.esys*wp.esys*x[i]*x[i]);
	      chi2+=x1;
	    }
	}
      else
	{
	  chi2=0;
	  for(i=1;i<=wp.npca;++i)
	    {
	      par_chi=0;
	      for(j=1;j<=wp.np;++j)
		par_chi+=wp.covar[j][i]*(x[j]-wp.x[j])/wp.e[j];
	      chi2+=par_chi*par_chi/wp.eigen[i];
	    }
	}
    }

      /* 
       * From Peder Norberg's instructions for use of PCA:

       do i=1,npca
       par_chi=0.d0
       do j=1,npoints
       par_chi=par_chi+pca(j,i)*(xi_teo(j)-xi_data(j))/err_data(j)
       enddo
       chi2=chi2+(par_chi**2)/ev(i)           ! (****)
       enddo

      */

  /* Add in the error on the galaxy density
   */
  if(HOD.free[0])
    chi2+=chi2ngal;

  t1 = clock();
  t0 = difftime(t1,t0)/CLOCKS_PER_SEC;
  niter++;
  if(!ThisTask){
    printf("ITER %7d %e ",niter,chi2);
    for(i=1;i<=ncf_hod;++i)printf("%e ",a[i]);
    printf(" %.2f\n",t0);
    fflush(stdout);
    if(HOD.free[0])
      printf("NGAL %e %e %e\n",chi2ngal,
	     GALAXY_DENSITY,wp.ngal);
  }
  chi2_prev=chi2;
  chi2_array[ichi]=chi2;
  ichi++;
  if(ichi==10)ichi=0;


  fflush(stdout);
  return(chi2);
}

void initial_wp_values(double *a, double **pp, double *yy)
{
  static int flag=1;
  int i,j;
  double d[100];


  if(flag) {
    i=0;j=0;
    if(HOD.free[++i])a[++j]=HOD.M_min;
    if(HOD.free[++i])a[++j]=HOD.M1;
    if(HOD.free[++i])a[++j]=HOD.alpha;
    if(HOD.free[++i])a[++j]=HOD.M_cut;
    if(HOD.free[++i])a[++j]=HOD.sigma_logM;
    if(HOD.free[++i])a[++j]=CVIR_FAC;
    if(HOD.pdfc>=7){
      if(HOD.free[++i])a[++j]=HOD.M_cen_max; }
    else {
      if(HOD.free[++i])a[++j]=HOD.MaxCen; }
    if(HOD.free[++i])a[++j]=HOD.M_sat_break;
    if(HOD.free[++i])a[++j]=HOD.alpha1;

    if(XCORR){
      i=0;
      if(HOD.free[++i])a[++j]=HOD2.M_min;
      if(HOD.free[++i])a[++j]=HOD2.M1;
      if(HOD.free[++i])a[++j]=HOD2.alpha;
      if(HOD.free[++i])a[++j]=HOD2.M_cut;
      if(HOD.free[++i])a[++j]=HOD2.sigma_logM;
      if(HOD.free[++i])a[++j]=CVIR_FAC;
      if(HOD2.pdfc>=7){
	if(HOD.free[++i])a[++j]=HOD2.M_cen_max; }
      else {
	if(HOD.free[++i])a[++j]=HOD2.MaxCen; }
    }
    printf("INITIAL VALUES: ");
    for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);   
    printf("\n");
  }
  //flag++;

  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=wp.ncf;++i)
    d[i]=a[i]*0.25/flag;


  if(POWELL)
    {
      for(i=1;i<=wp.ncf;++i)
	{
	  for(j=1;j<=wp.ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(j=1;j<=wp.ncf;++j)
	pp[1][j]=a[j];
      yy[1]=chi2_wp(a);
      
      for(i=1;i<=wp.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	  yy[i+1]=chi2_wp(a);	  
	  for(j=1;j<=wp.ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[wp.ncf]-=d[wp.ncf];
    }
}


/* This routine reads in the wp data and covariance matrix from
 * the filenames specified.
 * 
 * FILE FORMATS:
 * 
 * - fname_wp    ->  r xi e_xi
 * - fname_covar ->  (i=1,np)(j=1,np) read(covar[i][j])
 *
 */
void wp_input()
{
  float x1,x2,x3;
  FILE *fp;
  int i,j,n;
  char a[1000];

  if(!(fp=fopen(wp.fname_wp,"r")))
    {
      fprintf(stdout,"ERROR opening [%s]\n",wp.fname_wp);
      endrun("error in wp_input");
    }
  wp.np=filesize(fp);

  /* [wp.format==2] means that there are two header lines at 
   * the top of the file.
   */
  if(wp.format==2)
    {
      wp.np-=2;
      fgets(a,1000,fp);
      fgets(a,1000,fp);
    }

  /* [wp.format==3] means that there is one header lines at 
   * the top of the file.
   */
  if(wp.format==3)
    {
      wp.np-=1;
      fgets(a,1000,fp);
    }

  wp.r=dvector(1,wp.np);
  wp.x=dvector(1,wp.np);
  wp.e=dvector(1,wp.np);
  if(PCA)
    {
      wp.eigen=dvector(1,wp.np);
      wp.covar=dmatrix(1,wp.np,1,wp.np);
    }

  /* Read in the projected correlation function data.
   * Standard format [wp.format==1] is linear r, linear wp, linear err.
   * [wp.format==2] is log10 r, log10 wp, linear err.
   * NB! Peder's format is to list the inner edge of the bin, so we're adding 0.1 to each number.
   */
  for(i=1;i<=wp.np;++i)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      wp.r[i]=x1;
      wp.x[i]=x2;
      wp.e[i]=x3;
      if(wp.format==2){
	wp.r[i] = pow(10.0,wp.r[i]+0.1);
	wp.x[i] = pow(10.0,wp.x[i])*wp.r[i];
	wp.e[i] = wp.e[i]*wp.r[i];
      }
      if(wp.format==3){
	fscanf(fp,"%f",&x3);
	wp.e[i]=x3;
	wp.r[i] = pow(10.0,wp.r[i]+0.1);
	// wp.x[i] = wp.x[i]*wp.r[i];
	// wp.e[i] = wp.e[i]*wp.r[i];
      }
      if(wp.format==3 && PCA) {
	fscanf(fp,"%lf",&wp.eigen[i]);
	for(j=1;j<=wp.np;++j)
	  fscanf(fp,"%lf",&wp.covar[j][i]);
	if(wp.npca==0)
	  wp.npca = wp.np;
      }
	//      if(wp.format==3 && !PCA)
      fgets(a,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading %d lines from [%s]\n",wp.np,wp.fname_wp);

  if(!COVAR || PCA)
    return;
  /*
  if(wp.format==1)
    {
      Work.SDSS_bins=1;
      for(i=1;i<=40;++i)
	Work.rad[i] = i-0.5;
    }
  */
  if(!(fp=fopen(wp.fname_covar,"r")))
    {
      fprintf(stdout,"ERROR opening [%s]\n",wp.fname_covar);
      endrun("error in wp_input");
    }
  wp.covar=dmatrix(1,wp.np,1,wp.np);
  for(i=1;i<=wp.np;++i)
    for(j=1;j<=wp.np;++j)
      {
	fscanf(fp,"%lf",&(wp.covar[i][j]));
	/* printf("COVAR %d %d %e\n",i,j,wp.covar[i][j]); */
      }
  fclose(fp);
  if(!ThisTask)
    fprintf(stdout,"Done reading %d lines from [%s]\n",wp.np,wp.fname_covar);

}

