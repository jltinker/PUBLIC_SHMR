#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

/* Internal functions
 */
void initial_shmr_values(double *a, double **pp, double *yy);
double chi2_wrapper1(double *a);
int check_central_shmr_parameters(double *a);

/* External functions.
 */
double chi2_smf_shmr(void);
double chi2_wp_shmr(void);
int check_shmr_parameters(double *a);
void mcmc_clustering_input(void);

void shmr_minimization()
{
  static int icall = 1;
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m;
  FILE *fp;
  char aa[1000];

  CENTRALS_ONLY = 1;

  fprintf(stdout,"\n\nCHI2 MINIMIZATION OF W_P(R_P) DATA..........\n");
  fprintf(stdout,    "--------------------------------------------\n\n");

  if(POWELL)
    FTOL=1.0E-5;
  else
    FTOL=1.0E-5;

  n = SHMR_PARAMS + VARIABLE_EXCLUSION;
  if(CENTRALS_ONLY)
    {
      n = 5;
    }
  fprintf(stdout,"shmr_min> [%d] free parameters.\n",n);
  wpl.ncf = n;
  MCMC_OUTPUT = OUTPUT;

  mcmc_clustering_input();

  wpl.ncf=n;
  a=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  initial_shmr_values(a,pp,yy);


  if(POWELL) 
    {
      if(OUTPUT)printf("wp_min> starting powell.\n");
      powell(a,pp,n,FTOL,&niter,&chi2min,chi2_wrapper1);
      chi2min = chi2_wrapper1(a);
    }
  else
    {
      if(OUTPUT)printf("wp_min> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_wrapper1,&niter);
      for(i=1;i<=n;++i)a[i]=pp[1][i];
      chi2min = chi2_wrapper1(a);
    }	

  printf("POWELL %d %e ",icall++,chi2min);
  for(i=1;i<=n;++i)printf("%e ",a[i]);
  printf("\n");

  /* Output the fit 
   */
  sprintf(aa,"%s.fit",Task.root_filename);
  fp=fopen(aa,"w");
  fprintf(fp,"%e %e ",chi2min,HOD.M_min);
  for(i=1;i<=n;++i)fprintf(fp,"%e ",a[i]);
  fprintf(fp," %f\n",GALAXY_BIAS);
  fclose(fp);

  free_dvector(a,1,n);
  if(POWELL)
    free_dmatrix(pp,1,n,1,n);
  else
    free_dmatrix(pp,1,n+1,1,n);
  free_dvector(yy,1,n+1);

}

void initial_shmr_values(double *a, double **pp, double *yy)
{
  static int flag=1;
  int i,j;
  double d[100];

  muh(wpl.ncf);
  /* use initial inputs as startup values
   */
  for(i=1;i<=wpl.ncf;++i)
    a[i] = wpl.a[i];

  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=wpl.ncf;++i)
    d[i]=a[i]*0.1;

  /* smaller for first two variables.
   */
  d[1] /= 100.0;
  d[2] /= 100.0;

  /* If variable alpha, make d large
   */
  if(VARIABLE_ALPHA)
    {
      d[wpl.ncf-1] = -1.5*a[wpl.ncf-1];
      d[wpl.ncf] = -1.5*a[wpl.ncf];
    }
  

  if(POWELL)
    {
      for(i=1;i<=wpl.ncf;++i)
	{
	  for(j=1;j<=wpl.ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(j=1;j<=wpl.ncf;++j)
	pp[1][j]=a[j];
      yy[1]=chi2_wrapper1(a);
    
      for(i=1;i<=wpl.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	  yy[i+1]=chi2_wrapper1(a);	  
	  for(j=1;j<=wpl.ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[wpl.ncf]-=d[wpl.ncf];
    }
  if(OUTPUT)
    fprintf(stdout,"shmr_min> done with init\n");
}

double chi2_wrapper1(double *a)
{
  static int niter=0;
  int i;
  double chi2a, chi2b, chi2, t0, t1;
  for(i=1;i<=wpl.ncf;++i)
    wpl.a[i] = a[i];

  if(check_central_shmr_parameters(a)) { 
    muh(i=check_central_shmr_parameters(a)); 
    fmuh(a[abs(i)]);
    return 1.0E+7; }

  t0 = second();
  chi2a = chi2_wp_shmr();
  chi2b = chi2_smf_shmr();
  chi2 = chi2a+chi2b;
  t1 = second();
  
  fprintf(stdout,"ITER %d %e %e %e ",++niter,chi2a,chi2b,chi2);
  for(i=1;i<=wpl.ncf;++i)
    fprintf(stdout,"%e ",a[i]);
  fprintf(stdout,"%.2f\n",timediff(t0,t1));
  fflush(stdout);
  return chi2;
}

/* check and make sure the parameters selected are within the 
 * acceptable parameter space.
 */
int check_central_shmr_parameters(double *a)
{
  int check_flag, ibuf;
  float m1;

  wpl.reset_fred = 1;

  if(a[1]<=12.0)return 1;      // Mhalo_norm
  if(a[1]>=13.8)return -1;
    
  if(a[2]<=10.3)return 2;      // Mstellar norm
  if(a[2]>=12.0)return -2;
 
  if(a[3]<=0.2)return 3;       // Beta 
  if(a[3]>=1.8)return -3;

  if(a[4]<=0.01)return 4;       // delta (If <0, affects the threshold HOD, strange things happen .. don't go below 0.3)
  if(a[4]>=1.4)return -4;
  
  if(a[5]<=-0.5)return 5;      // gamma (can be <0)
  if(a[5]>=5.0)return -5;

  /*
  if(a[6]<=0.01)return 6;  // sigma_log_sm
  if(a[6]>=0.5)return -6;

  a[6] = HOD.sigma_logM;
  */
  return 0;
}
