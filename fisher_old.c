#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>

#include "header.h"

/* calculating expected errors on HOD+cosmological parameters for the NSF grant.
 */
void covar_pca_xx(double **covar, double **evect, double *eval, int ndata);


/* local functions
 */
double chi2rsd(double *a);
double chi2rsd_pca(double *a);
double chi2rsd_file(double *a);
double chi2wp(double *a);
double wp_input2(void);
double qp_input2(void);
double qp_input2a(void);
double qp_input3(void);

int FIT_RSD = 0;
int ICUT_RSD = 1000;
int ICUT_RSD_CUBE = 1000;
int FIT_WP = 1;
int ISKIP_RSD = 0;
int ITASK = 0;
int USE_FILES = 0;
int OUTPUT_DATAFILE = 0;
int FIT_QUADRUPOLE = 0;

double omz0, s8z0, sigv0;

void fisher()
{
  int i,j, nparams=4,n=4,subfish[100],i1,j1,ntemp,i2;
  double **tmp, **tmp2, **fisher, delta=1E-4, *a, *a0, *chi2xp, *chi2xm, **fcovar,
    **chi2xpp, **chi2xmm, *deltax,**ftemp, **chi2xpm, **chi2xmp;
  double x0,axis2,axis1,theta, xi_mono, xi_quad, x2,x1, gbias;
  FILE *fp;

// Define all the used matrices
  gsl_matrix * mxx = gsl_matrix_alloc (n, n);
  gsl_matrix * inverse = gsl_matrix_alloc (n, n);
  gsl_permutation * perm = gsl_permutation_alloc (n);
  int s;

  /*
  ITASK = -1;
  if(ARGC>3)
    ITASK = atoi(ARGV[3]);
  printf("ITASK %d\n",ITASK);
  if(ITASK==-1)USE_FILES = 1;
  if(ARGC>4)
    ISKIP_RSD = atoi(ARGV[4]);
  */

  nparams = 11;
  fisher = dmatrix(1,nparams,1,nparams);
  fcovar = dmatrix(1,nparams,1,nparams);
  a = dvector(1,nparams);
  a0 = dvector(1,nparams);
  deltax = dvector(1,nparams);
  chi2xp = dvector(1,nparams);
  chi2xm = dvector(1,nparams);
  chi2xpp = dmatrix(1,nparams,1,nparams);
  chi2xmm = dmatrix(1,nparams,1,nparams);
  chi2xmp = dmatrix(1,nparams,1,nparams);
  chi2xpm = dmatrix(1,nparams,1,nparams);


  // temp
  //HUBBLE = 0.7*(1+delta);

  //a0[1] = log10(HOD.M1);
  a0[1] = (HOD.M1);
  a0[2] = HOD.alpha;
  //a0[3] = log10(HOD.M_cut); //NB! there was once 100 here... putting it in hod.bat
  a0[3] = (HOD.M_cut); //NB! there was once 100 here... putting it in hod.bat
  a0[4] = HOD.sigma_logM;
  a0[5] = CVIR_FAC;
  a0[6] = HUBBLE;
  a0[7] = OMEGA_B;
  a0[8] = OMEGA_M;
  a0[9] = SIGMA_8;
  a0[10] = VBIAS;
  a0[11] = GAMMA;
  OMEGA_Z = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));

  //delta = 0.000001;
  delta = 1.0E-4;
  //delta = 1.0/atof(ARGV[2]);
  deltax[1] = delta;
  deltax[2] = delta;
  deltax[3] = delta;
  deltax[4] = delta;
  deltax[5] = delta;
  deltax[6] = delta;
  deltax[7] = delta;
  deltax[8] = delta;
  deltax[9] = delta;
  deltax[10] = delta;
  deltax[11] = delta;

  RESET_FLAG_1H++;
  RESET_FLAG_2H++;
  HOD.M_min = 0;
  set_HOD_params();
  //chi2wp(a0);
  gbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  omz0 = a0[8]*pow(1+REDSHIFT,3.0)/(a0[8]*pow(1+REDSHIFT,3.0)+(1-a0[8]));
  s8z0 = SIGMA_8;
  sigv0 = SIGV;
  BETA = pow(OMEGA_Z,GAMMA)/gbias;
  printf("BETA %f bias= %f\n",BETA,gbias);
  SIGV = sigv0*pow(OMEGA_Z/omz0,GAMMA)*(SIGMA_8/s8z0);
  //qp_input2();
  //qp_input2a();
  //qp_input3();
  

  //lets set up the multipoles too
  if(FIT_RSD){
    RESET_ZSPACE=1;
    //xi_multipoles(1.0,&xi_mono, &xi_quad);
  }

  /*
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;
  HOD.M_min = 0;
  set_HOD_params();
  */

  // assume that the input HOD is the proper HOD for the CMASS sample.
  // read in the "data" and the covariance matrix
  wp_input();

  if(!USE_FILES) {
  // replace the wp values with the ones from the analytic HOD
  for(i=1;i<=wp.np;++i)
    wp.x[i] = projected_xi(wp.r[i]);
  muh(wp.np);
  // replace the xi0 values with the ones from the analytic HOD
  if(FIT_RSD) {
  for(i=1;i<=rsdm.np;++i) {
    xi_multipoles(rsdm.r[i],&x0,&x2);
    //printf("%f %e %e %e %e\n",rsdm.r[i],rsdm.x[i],x0,
    //	   one_halo_real_space(rsdm.r[i])+two_halo_real_space(rsdm.r[i]),x2); 
    rsdm.x[i] = x0; rsdq.x[i] = x2; }
  }
  }

  if(OUTPUT_DATAFILE)
    {
      // lets output to a file for future use
      fp = fopen("wp_mockmean.dat","w");
      for(i=1;i<=wp.np;++i)
	fprintf(fp,"%e %e %e\n",wp.r[i],wp.x[i],sqrt(wp.covar[i][i]));
      fclose(fp);
      fp = fopen("wp_mockmean.covar","w");
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  fprintf(fp,"%e\n",wp.covar[i][j]);
      fclose(fp);
      /*
      fp = fopen("monopole_mockmean.dat","w");
      for(i=1;i<=rsdm.np;++i)
	fprintf(fp,"%e %e %e\n",rsdm.r[i],rsdm.x[i],rsdm.e[i]);
      fclose(fp);
      fp = fopen("monopole_mockmean.covar","w");
      for(i=1;i<=rsdm.np;++i)
	for(j=1;j<=rsdm.np;++j)
	  fprintf(fp,"%e\n",rsdm.covar[i][j]);
      fclose(fp);
      */
      fp = fopen("Q_mockmean.dat","w");
      for(i=1;i<=rsdq.np;++i)
	fprintf(fp,"%e %e %e\n",rsdq.r[i],rsdq.x[i],rsdq.e[i]);
      fclose(fp);
      fp = fopen("Q_mockmean.covar","w");
      for(i=1;i<=rsdq.np;++i)
	for(j=1;j<=rsdq.np;++j)
	  fprintf(fp,"%e\n",rsdq.covar[i][j]);
      fclose(fp);
      exit(0);
    }


  for(i=1;i<=wp.np;++i)
    {
      wp.e[i] = sqrt(wp.covar[i][i]);
    }
  //invert the matrix for wp
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
  
    //invert the matrix for rsdm
  if(FIT_RSD)
    {
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
    }

  if(FIT_QUADRUPOLE)
    {
      tmp=dmatrix(1,rsdq.np,1,1);
      tmp2=dmatrix(1,rsdq.np,1,rsdq.np);
      for(i=1;i<=rsdq.np;++i)
	for(j=1;j<=rsdq.np;++j)
	  tmp2[i][j]=rsdq.covar[i][j];
      gaussj(tmp2,rsdq.np,tmp,1);
      for(i=1;i<=rsdq.np;++i)
	for(j=1;j<=rsdq.np;++j)
	  rsdq.covar[i][j]=tmp2[i][j];
      free_dmatrix(tmp,1,rsdm.np,1,1);
      free_dmatrix(tmp2,1,rsdm.np,1,rsdm.np);
    }

  // swith to pca analysis for the chi2
  //covar_pca_xx(rsdm.covar, rsdm.evect, rsdm.eval, rsdm.np);
  //for(i=1;i<=rsdm.np;++i)
  // printf("COVAR %d %e\n",i,rsdm.eval[i]);

  // get the central chi2 (shouldn't this be 0 by definition?)
  x0 = chi2wp(a0) + chi2rsd(a0);
  x0 = 0;

  for(i=1;i<=nparams;++i)
    a[i] = a0[i];

  nparams = 4;
  
  // get the ii chi2 values
  for(i=1;i<=nparams;++i)
    {
      a[i] = a0[i]*(1-deltax[i]);
      printf("diff1 %d %e\n",i,a0[i]-a[i]);
      x1 = chi2wp(a);
      x2 = chi2rsd(a);
      chi2xp[i] = x1+x2;
      a[i] = a0[i]*(1-deltax[i]);
      printf("diff1 %d %e\n",i,a0[i]-a[i]);
      x1 = chi2wp(a);
      x2 = chi2rsd(a);
      chi2xm[i] = x1+x2;
      //chi2xm[i] = chi2xp[i];
      if(chi2xp[i]<x0)printf("ERRORp: %d %e %e\n",i,chi2xp[i],x0);
      if(chi2xm[i]<x0)printf("ERRORm: %d %e %e\n",i,chi2xm[i],x0);
      printf("CHI2XX %d %e %e %e %e %e %e\n",i,delta,chi2xp[i],chi2xm[i], 
	     (chi2xp[i] - 2*x0 + chi2xm[i])/(deltax[i]*deltax[i]*a0[i]*a0[i]),a[i],a0[i]);
      a[i] = a0[i];
    }

  //x1 = chi2wp(a);
  //x2 = chi2rsd(a);
  //fprintf(stdout,"ERROR CHECK1: %e %e\n",x1,x2);


  // get the off-diagonal chi2 values
  for(i=1;i<=nparams;++i)
    for(j=i+1;j<=nparams;++j)
      {
	a[i] = a0[i]*(1-deltax[i]);
	a[j] = a0[j]*(1-deltax[j]);
	x1 = chi2wp(a);
	x2 = chi2rsd(a);
	chi2xpp[i][j] = x1+x2;

	a[i] = a0[i]*(1-deltax[i]);
	a[j] = a0[j]*(1-deltax[j]);
	x1 = chi2wp(a);
	x2 = chi2rsd(a);
	chi2xmm[i][j] = x1+x2;
	//chi2xmm[i][j] = chi2xpp[i][j];

	a[i] = a0[i];
	a[j] = a0[j];
	continue;

	a[i] = a0[i]*(1+deltax[i]);
	a[j] = a0[j]*(1-deltax[j]);
	x1 = chi2wp(a);
	x2 = chi2rsd(a);
	chi2xpm[i][j] = x1+x2;

	a[i] = a0[i]*(1-deltax[i]);
	a[j] = a0[j]*(1+deltax[j]);
	x1 = chi2wp(a);
	x2 = chi2rsd(a);
	chi2xmp[i][j] = x1+x2;

      }

  //x1 = chi2wp(a);
  //x2 = chi2rsd(a);
  //fprintf(stdout,"ERROR CHECK2: %e %e\n",x1,x2);


  // now construct the fisher matrix diagonals
  for(i=1;i<=nparams;++i)
    fisher[i][i] = (chi2xp[i] - 2*x0 + chi2xm[i])/(deltax[i]*deltax[i]*a0[i]*a0[i]);

  // now construct the off-diagonal elements of the fisher matrix
  for(i=1;i<=nparams;++i)
    for(j=i+1;j<=nparams;++j)
      fisher[i][j] = (chi2xpp[i][j] - chi2xp[i] - chi2xp[j] + 2*x0 - chi2xm[i] - chi2xm[j] + chi2xmm[i][j])/
	(2*deltax[i]*deltax[j]*a0[i]*a0[j]);

  // now construct the off-diagonal elements of the fisher matrix (alternative method)
  for(i=1;i<=-nparams;++i)
    for(j=i+1;j<=nparams;++j)
      fisher[i][j] = (chi2xpp[i][j] - chi2xpm[i][j] - chi2xmp[i][j] + chi2xmm[i][j])/
	(4*deltax[i]*deltax[j]*a0[i]*a0[j]);

  for(i=1;i<=nparams;++i)
    for(j=i+1;j<=nparams;++j)
      fisher[j][i] = fisher[i][j];

  // put factor 1/2 in there
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      fisher[i][j]*=0.5;

  // get the subset of parameters we want to use.
  subfish[1] = 1;
  subfish[2] = 1;
  subfish[3] = 1;
  subfish[4] = 1;
  subfish[5] = 0;
  subfish[6] = 0;
  subfish[7] = 0;
  subfish[8] = 0;
  subfish[9] = 0;
  subfish[10] = 0;
  subfish[11] = 0;

  // lets get the sub-fisher-matrxi
  ntemp=0;
  for(i=1;i<=nparams;++i)
    if(subfish[i])ntemp++;

  ftemp = dmatrix(1,ntemp,1,ntemp);

  i1 = 0;
  for(i=1;i<=nparams;++i)
    {
      if(!subfish[i])continue;
      i1++;
      j1 = 0;
      a[i1] = a0[i];
      for(j=1;j<=nparams;++j)
	{
	  if(!subfish[j])continue;
	  j1++;
	  ftemp[i1][j1] = fisher[i][j];
	}
    }
  free_dmatrix(fisher,1,nparams,1,nparams);
  nparams = ntemp;
  fisher = dmatrix(1,ntemp,1,ntemp);
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      fisher[i][j] = ftemp[i][j];
  for(i=1;i<=nparams;++i)
    a0[i] = a[i];
			      
  // let's print out the fisher matrix
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      printf("FISHER %d %d %e\n",i,j,fisher[i][j]);

  // invert the fisher matrix for covariance matrix
  tmp=dmatrix(1,nparams,1,1);
  tmp2=dmatrix(1,nparams,1,nparams);
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      tmp2[i][j]=fisher[i][j];
  gaussj(tmp2,nparams,tmp,1);
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      fcovar[i][j]=tmp2[i][j];
  free_dmatrix(tmp,1,nparams,1,1);
  free_dmatrix(tmp2,1,nparams,1,nparams);
  
  
  
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      gsl_matrix_set(mxx,i-1,j-1,fisher[i][j]);
  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (mxx, perm, &s);
  // Invert the matrix m
  gsl_linalg_LU_invert (mxx, perm, inverse);
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      printf("GSL %d %d %e %e\n",i,j,fcovar[i][j], gsl_matrix_get(inverse,i-1,j-1));
  //      fcovar[i][j] = gsl_matrix_get(inverse,i-1,j-1);
  

  // test: matrix multiplication
  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      {
	x1 = 0;
	for(n=1;n<=nparams;++n)
	  x1+=fcovar[n][i]*fisher[j][n];
	//x1+=fcovar[i][n]*fisher[n][j];
	printf("MULT %d %d %e\n",i,j,x1);
      }

  for(i=1;i<=nparams;++i)
    for(j=1;j<=nparams;++j)
      printf("%d %d %e %e %e\n",i,j,fisher[i][j],fcovar[i][j],fcovar[i][j]/sqrt(fabs(fcovar[i][i]*fcovar[j][j])));
  for(i=1;i<=-nparams;++i)
    for(j=1;j<=nparams;++j)
      fcovar[i][j] = fabs(fcovar[i][j]);

  // lets print out the diagonal elements:
  for(i=1;i<=nparams;++i)
    printf("a[%d]=%e, c= %e, frac= %e  %e\n",i,a0[i],sqrt(fabs(fcovar[i][i])),
	   sqrt(fabs(fcovar[i][i]))/a0[i],fcovar[i][i]);///abs(fcovar[i][i]));

  //lets print out the elements of the ellipses
  for(i=1;i<=nparams;++i)
    for(j=i+1;j<=nparams;++j)
      {
	axis1 = 0.5*(fcovar[i][i] + fcovar[j][j]) + sqrt((fcovar[i][i]-fcovar[j][j])*(fcovar[i][i]-fcovar[j][j])*0.25 + 
						     fcovar[i][j]*fcovar[i][j]);
	axis2 = 0.5*(fcovar[i][i] + fcovar[j][j]) - sqrt((fcovar[i][i]-fcovar[j][j])*(fcovar[i][i]-fcovar[j][j])*0.25 + 
						     fcovar[i][j]*fcovar[i][j]);
	theta = 0.5*atan(2*fcovar[i][j]/(fcovar[i][i]-fcovar[j][j]))*180/PI;
	axis2 = fabs(axis2);
	printf("%d %d %e %e %e %e %e %e\n",
	       i,j,sqrt(axis1),sqrt(axis2),theta,sqrt(fcovar[i][i]),sqrt(fcovar[j][j]),fcovar[i][j]);
      }
  exit(0);
}

double chi2wp(double *a)
{
  static int niter=0;
  int i,j;
  double x1, chi2 = 0, x[100], gbias;
  FILE *fp;

  if(USE_FILES)return 0;
  
  /*  
  if((niter+1) % 12 != ITASK)
    {
      niter++;
      return 0;
    }
  */
  printf("Calculating for %d\n",niter+1);

  // put the parameters in place
  //HOD.M1 = pow(10.0,a[1]);
  HOD.M1 = a[1];
  HOD.alpha = a[2];
  //HOD.M_cut = pow(10.0,a[3]);
  HOD.M_cut = a[3];
  HOD.sigma_logM = a[4];
  CVIR_FAC = a[5];
  //VBIAS = a[6];
  //VBIAS_C = a[7];
  HUBBLE = a[6];
  OMEGA_B = a[7];

  // and cosmology;
  OMEGA_M = a[8];
  OMEGA_Z = OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  SIGMA_8 = a[9];
  VBIAS = a[10];
  GAMMA = a[11];

  GALAXY_DENSITY = 0.0003;
  HOD.M_min = 0;
  set_HOD_params();
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;


  gbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  BETA = pow(OMEGA_Z,GAMMA)/gbias;
  SIGV = sigv0*pow(OMEGA_Z/omz0,GAMMA)*(SIGMA_8/s8z0);
  printf("BETA %f bias= %f\n",BETA,gbias);
  //RESET_COSMOLOGY++;


  //printf("%e %e %e %e %e\n",HOD.M_min,a[1],a[2],a[3],a[4]);
  if(!FIT_WP)return 0;

  for(i=1;i<=wp.np;++i)
    x[i] = projected_xi(wp.r[i]);



  for(i=1;i<=wp.np;++i)
    printf("WP %d %e %e %e %e\n",niter+1,wp.r[i],wp.x[i],wp.e[i],x[i]);

  //let's try diagonals
  for(i=1;i<=wp.np;++i)
    chi2 += (x[i]-wp.x[i])*(x[i]-wp.x[i])/wp.e[i]/wp.e[i];


  for(i=1;i<=-wp.np;++i)
    for(j=1;j<=wp.np;++j)
      {
	//if(i==j)printf("%d %e %e %e\n",i,wp.r[i],wp.x[i],x[i]);
	x1=(x[i]-wp.x[i])*(x[j]-wp.x[j])*wp.covar[i][j];
	chi2+=x1;
      }
  printf("chi2wp> %d %e\n",++niter,chi2);
  return chi2;
}

double chi2rsd(double *a)
{
  static int niter=0;
  int i,j;
  double x1, x2,chi2 = 0, x[100];
  FILE *fp;

  if(!FIT_RSD) return 0;
  if(USE_FILES)return chi2rsd_file(a);
  //if(!USE_FILES)return chi2rsd_pca(a);
  //return chi2rsd_file(a);
  //if(rsdm.npca)return chi2rsd_pca(a);
  //if(!FIT_RSD)return 0;
  
  RESET_PVZ = 1;
  RESET_ZSPACE = 1;

  for(i=1;i<=rsdm.np;++i)
    xi_multipoles(rsdm.r[i],&x[i],&x2);

  // lets output to a file for future use

  for(i=1;i<=rsdm.np;++i)
    for(j=1;j<=rsdm.np;++j)
      {
	if(i==j)printf("rsd%d %d %e %e %e %e\n",niter+1,i,rsdm.r[i],rsdm.x[i],x[i],rsdm.e[i]);
	x1=(x[i]-rsdm.x[i])*(x[j]-rsdm.x[j])*rsdm.covar[i][j];
	chi2+=x1;
      }
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;
}

double chi2rsdq(double *a)
{
  static int niter=0;
  int i,j;
  double x1, x2,chi2 = 0, x[100];
  FILE *fp;

  RESET_PVZ = 1;
  RESET_ZSPACE = 1;

  for(i=1;i<=rsdq.np;++i)
    xi_multipoles(rsdq.r[i],&x1,&x[i]);

  // if niter==0, then set this as fiducial model
  if(niter==0)
    {
      for(i=1;i<=rsdq.np;++i)
	rsdq.x[i] = x[i];
    }

  // lets output to a file for future use
  for(i=1;i<=rsdq.np;++i)
    for(j=1;j<=rsdq.np;++j)
      {
	if(i==j)printf("rsd%d %d %e %e %e %e\n",niter+1,i,rsdq.r[i],rsdq.x[i],x[i],rsdq.e[i]);
	x1=(x[i]-rsdq.x[i])*(x[j]-rsdq.x[j])*rsdq.covar[i][j];
	chi2+=x1;
      }
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;
}

double chi2rsd_file(double *a)
{
  static int niter = 0;
  int i,j,k;
  FILE *fp;
  char fname[1000];
  float x1;
  double chi2 =0, x[100], x2,xq[100];

  muh(niter);

  if(!FIT_RSD)return 0;
  
  // get the fiducial model  
  fp = openfile("model.1");
  if(ISKIP_RSD)
    for(i=1;i<=ISKIP_RSD;++i)
      fscanf(fp,"%s %d %e %e %e %lf %lf",fname,&k,&x1,&x1,&x1,&x2,&x2);      
  for(i=1;i<=rsdm.np;++i)
    fscanf(fp,"%s %d %e %e %e %lf %lf",fname,&k,&x1,&x1,&x1,&rsdm.x[i],&rsdq.x[i]);
  fclose(fp);

  //get the current model
  sprintf(fname,"model.%d",niter+1);
  fp = openfile(fname);
  if(ISKIP_RSD)
    for(i=1;i<=ISKIP_RSD;++i)
      fscanf(fp,"%s %d %e %e %e %lf %lf",fname,&k,&x1,&x1,&x1,&x2,&x2);      
  for(i=1;i<=rsdm.np;++i)
    fscanf(fp,"%s %d %e %e %e %lf %lf",fname,&k,&x1,&x1,&x1,&x[i],&xq[i]);
  fclose(fp);

  if(FIT_QUADRUPOLE)
    {
      // get the chi2
      for(i=1;i<=rsdm.np;++i)
	for(j=1;j<=rsdm.np;++j)
	  {
	    if(i==j)printf("RSD%d %d %e %e %e %e\n",niter+1,i,rsdq.r[i],rsdq.x[i],xq[i],rsdq.e[i]);
	    x1=(xq[i]-rsdq.x[i])*(xq[j]-rsdq.x[j])*rsdq.covar[i][j];
	    chi2+=x1;
	  }
      //printf("chi2rsd> %d %e\n",++niter,chi2);
      //return chi2;
    }

  // get the chi2
  for(i=1;i<=rsdm.np;++i)
    for(j=1;j<=rsdm.np;++j)
      {
	if(i==j)printf("RSD%d %d %e %e %e\n",niter+1,i,rsdm.r[i],rsdm.x[i],x[i]);
	x1=(x[i]-rsdm.x[i])*(x[j]-rsdm.x[j])*rsdm.covar[i][j];
	chi2+=x1;
      }
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;

  // DIAG ERRORS!
  chi2 = 0;
  for(i=1;i<=rsdm.np;++i) 
    chi2 += (x[i]-rsdm.x[i])*(x[i]-rsdm.x[i])/(rsdm.e[i]*rsdm.e[i]); 
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;


}

double chi2rsd_pca(double *a)
{
  double chi1,chi2;
  int i,j,k,n;

  static int niter=0;
  double x1, x2, x[100], q[100];

  if((niter+1) % 12 != ITASK)
    {
      niter++;
      return 0;
    }
  printf("Calculating for %d\n",niter+1);

  if(!FIT_RSD)return 0;
  
  RESET_PVZ = 1;
  RESET_ZSPACE = 1;

  printf("CHIP%d ",niter+1);
  for(i=1;i<=9;++i)
    printf("%e ",a[i]);
  printf("\n");

  for(i=1;i<=rsdm.np;++i)
    xi_multipoles(rsdm.r[i],&x[i],&q[i]);
  for(i=1;i<=rsdm.np;++i)
    printf("RSD %d %e %e %e %.16e %.16e\n",niter+1,rsdm.r[i],rsdm.x[i],rsdm.e[i],x[i],q[i]);
 
  chi2 = 0;
  for(i=1;i<=rsdm.npca;++i)
    {
      chi1 = 0;
      for(j=1;j<=rsdm.np;++j)
	chi1 += rsdm.evect[j][i]*(x[j]-rsdm.x[j])/rsdm.e[j];
      chi2 += chi1*chi1/rsdm.eval[i];
    }
  
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;

 // DIAG ERRORS!
  chi2 = 0;
  for(i=1;i<=rsdm.np;++i) 
    chi2 += (x[i]-rsdm.x[i])*(x[i]-rsdm.x[i])/(rsdm.e[i]*rsdm.e[i]); 
  printf("chi2rsd> %d %e\n",++niter,chi2);
  return chi2;
 
}


double wp_input2()
{
  int i,j,k,n;
  char fname[1000], aa[1000];
  float x1,x2,x3,xx[100];
  FILE *fp;

  muh(0);

  // get the number of data points;
  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_WP/wp_%04d.dat",1);
  fp = openfile(fname);
  wp.np = filesize(fp);
  wp.covar = dmatrix(1,wp.np,1,wp.np);
  wp.r = dvector(1,wp.np);
  wp.x = dvector(1,wp.np);
  wp.e = dvector(1,wp.np);

  for(i=1;i<=wp.np;++i)
    for(j=1;j<=wp.np;++j)
      wp.covar[i][j] = 0;

  n = 1000;
  for(j=1;j<=wp.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      wp.x[j] = x3/n;
      wp.r[j] = x1; 
    }
  fclose(fp);

  for(i=2;i<=n;++i)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_WP/wp_%04d.dat",i);
      fp = openfile(fname);
      for(j=1;j<=wp.np;++j)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  fgets(aa,1000,fp);
	  wp.x[j] += x3/n;
	}
      fclose(fp);
      //muh(i);
    }

  for(k=1;k<=n;++k)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_WP/wp_%04d.dat",k);
      fp = openfile(fname);
      for(i=1;i<=wp.np;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&xx[i]);
	  fgets(aa,1000,fp);
	}
      fclose(fp);
	    
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	    wp.covar[i][j] += (xx[i]-wp.x[i])*(xx[j]-wp.x[j])/n;
    }

  for(i=1;i<=wp.np;++i)
    printf("WP %e %e %e\n",wp.r[i],wp.x[i],sqrt(wp.covar[i][i]));
  //exit(0);

     
}

double qp_input2()
{
  int i,j,k,n,nfiles=1000,itail;
  char fname[1000], aa[1000], froot[100], fname1[1000], command[1000];
  float x1,x2,x3,xx[100],xibar;
  FILE *fp;

  sprintf(froot,"ximulti");

  if(ISKIP_RSD)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
      fp = openfile(fname);
      rsdm.np = filesize(fp);
      fclose(fp);
      sprintf(froot,"ximulti_xx");
      itail = rsdm.np - ISKIP_RSD;
      for(i=1;i<=nfiles;++i)
	{
	  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/ximulti_%04d.dat",i);
	  sprintf(fname1,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
	  sprintf(command,"tail -%d %s > %s",itail,fname,fname1);
	  fprintf(stdout,"input_rsd> [%s]\n",command);
	  system(command);
	} 
    }

  // get the number of data points;
  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
  fp = openfile(fname);
  rsdm.np = filesize(fp);
  if(rsdm.np>ICUT_RSD)rsdm.np = ICUT_RSD;
  rsdm.covar = dmatrix(1,rsdm.np,1,rsdm.np);
  rsdm.r = dvector(1,rsdm.np);
  rsdm.x = dvector(1,rsdm.np);
  rsdm.e = dvector(1,rsdm.np);



  for(i=1;i<=rsdm.np;++i)
    for(j=1;j<=rsdm.np;++j)
      rsdm.covar[i][j] = 0;

  n = nfiles;
  for(j=1;j<=rsdm.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      rsdm.x[j] = x2/n;
      rsdm.r[j] = x1; 
    }
  fclose(fp);


  for(i=2;i<=n;++i)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
      fp = openfile(fname);
      for(j=1;j<=rsdm.np;++j)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  fgets(aa,1000,fp);
	  rsdm.x[j] += x2/n;
	}
      fclose(fp);
      //muh(i);
    }

  for(k=1;k<=n;++k)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,k);
      fp = openfile(fname);
      for(i=1;i<=rsdm.np;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&xx[i],&x2);
	  fgets(aa,1000,fp);
	}
      fclose(fp);
	    
      for(i=1;i<=rsdm.np;++i)
	for(j=1;j<=rsdm.np;++j)
	    rsdm.covar[i][j] += (xx[i]-rsdm.x[i])*(xx[j]-rsdm.x[j])/n;
    }
  for(i=1;i<=rsdm.np;++i)
    rsdm.e[i] = sqrt(rsdm.covar[i][i]);

  for(i=1;i<=-rsdm.np;++i)
    printf("INPUT %e %e %e\n",rsdm.r[i],rsdm.x[i],rsdm.e[i]);

  // lets do some PCA
  rsdm.eval = dvector(1,rsdm.np);
  rsdm.evect = dmatrix(1,rsdm.np,1,rsdm.np);
  rsdm.npca = 50;
      
}
double qp_input2a()
{
  int i,j,k,n,nfiles=1000,itail;
  char fname[1000], aa[1000], froot[100], fname1[1000], command[1000];
  float x1,x2,x3,xx[100],xibar,xb[100],rr;
  FILE *fp;

  sprintf(froot,"ximulti");

  if(ISKIP_RSD)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
      fp = openfile(fname);
      rsdq.np = filesize(fp);
      fclose(fp);
      sprintf(froot,"ximulti_xx");
      itail = rsdq.np - ISKIP_RSD;
      for(i=1;i<=nfiles;++i)
	{
	  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/ximulti_%04d.dat",i);
	  sprintf(fname1,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
	  sprintf(command,"tail -%d %s > %s",itail,fname,fname1);
	  fprintf(stdout,"input_rsd> [%s]\n",command);
	  system(command);
	} 
    }

  // get the number of data points;
  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
  fp = openfile(fname);
  rsdq.np = filesize(fp);
  if(rsdq.np>ICUT_RSD)rsdq.np = ICUT_RSD;
  rsdq.covar = dmatrix(1,rsdq.np,1,rsdq.np);
  rsdq.r = dvector(1,rsdq.np);
  rsdq.x = dvector(1,rsdq.np);
  rsdq.e = dvector(1,rsdq.np);



  for(i=1;i<=rsdq.np;++i)
    for(j=1;j<=rsdq.np;++j)
      rsdq.covar[i][j] = 0;

  n = nfiles;
  xibar = 0;
  for(j=1;j<=rsdq.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      xibar += x1*x1*x2;
      rr = j+ISKIP_RSD;
      rsdq.x[j] = x3/(x2-3./(rr*rr*rr)*xibar)/n;
      rsdq.r[j] = x1; 
      printf("BB %e %e %e %e\n",rsdq.r[j],x2,xibar,rsdq.x[j]*n);
    }
  fclose(fp);


  for(i=2;i<=n;++i)
    {
      xibar = 0;
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
      fp = openfile(fname);
      for(j=1;j<=rsdq.np;++j)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  fgets(aa,1000,fp);
	  xibar += x1*x1*x2;
	  rr = j+ISKIP_RSD;
	  rsdq.x[j] += x3/(x2-3./(rr*rr*rr)*xibar)/n;
	  rsdq.r[j] = x1; 
	}
      fclose(fp);
      //muh(i);
    }

  for(k=1;k<=n;++k)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,k);
      fp = openfile(fname);
      xibar = 0;
      for(i=1;i<=rsdq.np;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  xibar += x1*x1*x2;
	  rr = i+ISKIP_RSD;
	  xx[i] = x3/(x2-3./(rr*rr*rr)*xibar);
	  fgets(aa,1000,fp);
	  if(k==1)
	    printf("%e %e\n",x1,xx[i]);
	}
      fclose(fp);
	    
      for(i=1;i<=rsdq.np;++i)
	for(j=1;j<=rsdq.np;++j)
	    rsdq.covar[i][j] += (xx[i]-rsdq.x[i])*(xx[j]-rsdq.x[j])/n;
    }
  for(i=1;i<=rsdq.np;++i)
    rsdq.e[i] = sqrt(rsdq.covar[i][i]);

  for(i=1;i<=rsdq.np;++i)
    printf("INPUTq %e %e %e\n",rsdq.r[i],rsdq.x[i],rsdq.e[i]);

  // lets do some PCA
  rsdq.eval = dvector(1,rsdq.np);
  rsdq.evect = dmatrix(1,rsdq.np,1,rsdq.np);
  rsdq.npca = 50;
      
}



double qp_input3()
{
  int i,j,k,n,nfiles=200,itail,istart=100;
  char fname[1000], aa[1000], froot[100], fname1[1000], command[1000];
  float x1,x2,x3,xx[100],xibar;
  FILE *fp;

  sprintf(froot,"ximulti1");

  if(ISKIP_RSD)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/%s_%03d.dat",froot,istart);
      fp = openfile(fname);
      rsdm.npc = filesize(fp);
      fclose(fp);
      sprintf(froot,"ximulti1_xx");
      itail = rsdm.npc - ISKIP_RSD;
      for(i=istart;i<=istart+nfiles-1;++i)
	{
	  sprintf(fname,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/ximulti1_%03d.dat",i);
	  sprintf(fname1,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/%s_%03d.dat",froot,i);
	  sprintf(command,"tail -%d %s > %s",itail,fname,fname1);
	  fprintf(stdout,"input_rsd> [%s]\n",command);
	  system(command);
	} 
    }

  // get the number of data points;
  sprintf(fname,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/%s_%03d.dat",froot,istart);
  fp = openfile(fname);
  rsdm.npc = filesize(fp);
  if(rsdm.npc>ICUT_RSD_CUBE)rsdm.npc = ICUT_RSD_CUBE;
  rsdm.covarc = dmatrix(1,rsdm.npc,1,rsdm.npc);
  rsdm.rc = dvector(1,rsdm.npc);
  rsdm.xc = dvector(1,rsdm.npc);
  rsdm.ec = dvector(1,rsdm.npc);

  for(i=1;i<=rsdm.npc;++i)
    for(j=1;j<=rsdm.npc;++j)
      rsdm.covarc[i][j] = 0;

  n = nfiles;
  for(j=1;j<=rsdm.npc;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      rsdm.xc[j] = x2/n;
      rsdm.rc[j] = x1; 
    }
  fclose(fp);


  for(i=1+istart;i<=n+istart-1;++i)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/%s_%03d.dat",froot,i);
      fp = openfile(fname);
      for(j=1;j<=rsdm.npc;++j)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  fgets(aa,1000,fp);
	  rsdm.xc[j] += x2/n;
	}
      fclose(fp);
    }

  for(k=istart;k<=n+istart-1;++k)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/CUBIC_MOCKS/%s_%03d.dat",froot,k);
      fp = openfile(fname);
      for(i=1;i<=rsdm.npc;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&xx[i],&x2);
	  fgets(aa,1000,fp);
	}
      fclose(fp);
	    
      for(i=1;i<=rsdm.npc;++i)
	for(j=1;j<=rsdm.npc;++j)
	    rsdm.covarc[i][j] += (xx[i]-rsdm.xc[i])*(xx[j]-rsdm.xc[j])/n;
    }
  for(i=1;i<=rsdm.npc;++i)
    rsdm.ec[i] = sqrt(rsdm.covarc[i][i]);

  // add the diagonals
  for(i=1;i<=rsdm.npc;++i)
    rsdm.e[i] = sqrt(rsdm.covarc[i][i] + rsdm.covar[i][i]);


  for(i=1;i<=rsdm.np;++i)
    printf("INPUT %e %e %e %e %e\n",rsdm.r[i],rsdm.x[i],rsdm.e[i],rsdm.ec[i]/rsdm.x[i],sqrt(rsdm.covar[i][i])/rsdm.x[i]);
  exit(0);

  // combine the two covariance matrices
  for(i=1;i<=rsdm.npc;++i)
    for(j=1;j<=rsdm.npc;++j)
      rsdm.covar[i][j] += rsdm.covarc[i][j];


}

