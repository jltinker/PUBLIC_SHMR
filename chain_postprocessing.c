#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

/*------------------------------------------------------------------------
 * POST-PROCESSING OF A FULL CHAIN 
 * ----------------------------------------------------------------------*/

void chain_postprocess()
{
  int i,j,k,n,i1,nm,ibuf=11, FLAG;
  double *a, **red_sat_smf, **red_cen_smf, **sat_smf, **cen_smf,*xx, ** all_smf,
    **mhalo, **mhalo_red, mass, *mass_array, ngalr, ngalb;
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
  mlo = pow(10,11.4);
  if(REDSHIFT<0.2)
    mlo = 1.0E+10;
  mhi = pow(10.0,12.9);
  dlogm = log(mhi/mlo)/(nm-1);

  red_sat_smf = dmatrix(1,nm,1,n);
  red_cen_smf = dmatrix(1,nm,1,n);
  sat_smf = dmatrix(1,nm,1,n);
  cen_smf = dmatrix(1,nm,1,n);
  all_smf = dmatrix(1,nm,1,n);
  mhalo = dmatrix(1,nm,1,n);
  mhalo_red =  dmatrix(1,nm,1,n);

  // temp vector
  xx = dvector(1,n);

  // loop through all elements
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%d %d",&i1,&i1);
      for(j=1;j<=wpl.ncf;++j)
	fscanf(fp,"%lf",&a[j]);
      fgets(aa,1000,fp);

      // set up global
      for(j=1;j<=wpl.ncf;++j)
	wpl.a[j] = a[j];
      if(VARIABLE_EXCLUSION)
	EXCLUSION_RADIUS = a[wpl.ncf];

      fprintf(stderr,"post_process> %d\n",i);

      // now loop through all masses
      for(j=1;j<=nm;++j)
	{
	  mass = exp((j-1)*dlogm)*mlo;
	  mhalo[j][i] = ms_to_mhalo(mass,a);
	  
	  // get SMFs for blue gals
	  FLAG = 1;
	  wpl.reset_inversion = 1;
	  set_up_hod_for_shmr(exp(-dlogm/2)*mass, exp(dlogm/2)*mass, a);
	  ngalb = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
	  cen_smf[j][i] = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt)/dlogm*log(10)*pow(HUBBLE,3.0);
      	  sat_smf[j][i] = ngalb-cen_smf[j][i];
	  all_smf[j][i] = ngalb;
	}
    }

  sprintf(fname,"%s_stats.all_smf",Task.root_filename);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,all_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,all_smf[i][n/2],
	     all_smf[i][(int)(n*0.16)],all_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"%s_stats.cen_smf",Task.root_filename);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,cen_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,cen_smf[i][n/2],
	     cen_smf[i][(int)(n*0.16)],cen_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);

  sprintf(fname,"%s_stats.sat_smf",Task.root_filename);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      sort(n,sat_smf[i]);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,sat_smf[i][n/2],
	     sat_smf[i][(int)(n*0.16)],sat_smf[i][(int)(n*0.84)]);
    }
  fclose(fpout);


  sprintf(fname,"%s_stats.fsat",Task.root_filename);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = sat_smf[i][j]/(sat_smf[i][j]+cen_smf[i][j]);
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);


  sprintf(fname,"%s_stats.shmr",Task.root_filename);
  fpout = fopen(fname,"w");
  for(i=1;i<=nm;++i)
    {
      for(j=1;j<=n;++j)
	xx[j] = mhalo[i][j];
      sort(n,xx);
      fprintf(fpout,"%e %e %e %e\n",exp((i-1)*dlogm)*mlo,xx[n/2],
	     xx[(int)(n*0.16)],xx[(int)(n*0.84)]);
    }
  fclose(fpout);

}
