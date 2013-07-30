#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

void svdcmp(double **a, int m, int n, double w[], double **v);

void covar_pca(int ii)
{
  double **tmp, **tmp1, **evect, *eval, **tmp2, **diag, **tmp3;
  int n,i,j,k,nrot;
  double cumu = 0;

  n = wpl.ndata[ii];
  tmp = dmatrix(1,n,1,n);
  tmp2 = dmatrix(1,n,1,n);
  tmp3 = dmatrix(1,n,1,n);
  diag = dmatrix(1,n,1,n);
  tmp1 = dmatrix(1,n,1,1);
  evect = dmatrix(1,n,1,n);
  eval = dvector(1,n);

  wpl.evect[ii] = dmatrix(1,n,1,n);
  wpl.eval[ii] = dvector(1,n);
  wpl.npca = PCA;

  for(i=1;i<=n;++i)
    wpl.edata[ii][i] = sqrt(wpl.covar[ii][i][i]);


  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      tmp[i][j] = wpl.covar[ii][i][j]/wpl.edata[ii][i]/wpl.edata[ii][j];
  for(i=1;i<=-n;++i)
    for(j=1;j<=n;++j)
      tmp3[i][j] = tmp[i][j] = wpl.covar[ii][i][j];


  //jacobi(tmp,n,eval,evect,&nrot);
  /* replace jacobi with a singular-value decomposition
   */
  svdcmp(tmp,n,n,eval,evect);

  for(i=1;i<=n;++i)
    {
      wpl.eval[ii][i] = eval[i];
      for(j=1;j<=n;++j)
	wpl.evect[ii][i][j] = evect[i][j];
    }
  return;
}

double pca_chi2(double *xm, int ii)
{
  double chi1,chi2;
  int i,j,k,n;
  n = wpl.ndata[ii];

  chi2 = 0;
  for(i=1;i<=wpl.npca;++i)
    {
      chi1 = 0;
      for(j=1;j<=n;++j)
	chi1 += wpl.evect[ii][j][i]*(xm[j]-wpl.xdata[ii][j])/wpl.edata[ii][j];
      chi2 += chi1*chi1/wpl.eval[ii][i];
    }
  return chi2;
}
