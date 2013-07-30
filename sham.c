#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

double *mass, *decon, shambias[10], shambias2[10];
int nmass;

void convlv(double data[], unsigned long n, double respns[], unsigned long m,
	    int isign, double ans[]);
void bias_from_simulation_shmr();

void deconvolved_smf()
{
  static int icall=0;
  static double *smf, *rfunc;
  int i,j,k,n,m,npad,ibreak, flag = 1;
  double dm, mlo, mhi, sigma = 0.02, smf0;
  double a[10], xx;

  sigma = wpl.a[6];

  n = pow(2,25);
  nmass = n;
  mhi = 175.0;
  mlo = -200.0;
  dm = (mhi-mlo)/(n-1);
  // set up the response function
  m = (int)(4*sigma/dm*2);
  if(m%2==0) m++;
  fprintf(stderr,"M: %d %d\n",m,n);
  npad = m/2+1;
  fprintf(stderr,"NPAD: %d\n",npad);

  if(!icall)
    {
      icall++;
      mass = dvector(1,n);
      smf = dvector(1,n);
      rfunc = dvector(1,n);
      decon = dvector(1,2*n);
    }

  // old fit (fixing low-M slope to Main)
  a[1]=  1.406182e-01;
  a[2]=  1.203550e+01;
  a[3]=  -7.479410e-01;
  a[4]=  5.842199e-01;
  a[5]=  -8.706507e+00;

  // new fit (free all)
  
  a[1] = 2.430359e-02;
  a[2] = 1.146781e+01;
  a[3] = 7.280067e-01;
  a[4] = 1.172995e+00;
  a[5] = -6.012479e+00;
  
  for(i=1;i<=n;++i)
    {
      rfunc[i] = 0;
      mass[i] = mlo + (i-1)*dm;
      // use the analytic model fit
      xx = pow(10.0,(mass[i]-a[2]));
      smf[i] = a[1]*pow(xx,a[3]+1)*pow(1+pow(xx,a[4]),(a[5]-a[3])/a[4])*log(10);
      smf[i] += 0*1.0E-30;
      /*
      if(mass[i]>13 && flag) {
	smf0 = smf[i]/(a[1]*pow(xx,a[3]+1)*pow(1+pow(xx,a[4]),(-5-a[3])/a[4])*log(10));
	flag = 0; }
      if(mass[i]>13)
	smf[i] = smf0*a[1]*pow(xx,a[3]+1)*pow(1+pow(xx,a[4]),(-5-a[3])/a[4])*log(10);
      */
      if(isnan(smf[i]))smf[i] = 0.0;
      if(i>n-npad)smf[i] = 0;
    }

  for(i=1;i<=(m-1)/2+1;++i)
    rfunc[i] = 1/RT2PI/sigma*exp(-(mass[i]-mlo)*(mass[i]-mlo)/(2*sigma*sigma))*dm;
  for(j=1,i=m;i>(m-1)/2;--i,++j)
    rfunc[i] = rfunc[j];

  for(i=1;i<=-m;++i)
    printf("RFUNC %d %e\n",i,rfunc[i]);

  convlv(smf,n,rfunc,m,-1,decon);
  fprintf(stderr,"Done with deconvolution...\n");
  
  for(i=1;i<=n;++i)
    if(mass[i]>11.6)break;
  ibreak = i;

  for(i=npad;i<=n-npad;++i)
    if(i%64==0)
      printf("DECON %e %e %e %e %e\n",mass[i],smf[i],decon[i],smf[i]*decon[ibreak]/smf[ibreak],
	     smf[i+m/16]);
  fflush(stdout);
}

void bias_from_simulation()
{
  int nhalo, i,j,k, isat, nbins=8, ibin, jstart, *id;
  double xbias[10], hcnt[10], mbin[10];
  double m200, m, mpeak, *mhost, vpeak, mgal, mstar, *mx, volume, ncumu, density, dlogm, *x, *y, *z;
  FILE *fp;
  char aa[1000];
  unsigned long IDUM = 555;

  volume = (1000.0*1000.0*1000.0); // in h=1 units

  fp = openfile("halosubs_multidark_z0.53_D200.dat");
  nhalo = filesize(fp);

  mx = dvector(1,nhalo);
  mhost = dvector(1,nhalo);
  id = ivector(1,nhalo);
  x = dvector(1,nhalo);
  y = dvector(1,nhalo);
  z = dvector(1,nhalo);

  for(i=1;i<=nhalo;++i)
    {
      fscanf(fp,"%lf %lf %lf %lf %d %lf %lf %lf",
	     &m200, &mhost[i], &mpeak, &vpeak, &isat, &x[i], &y[i], &z[i]);
      fgets(aa,1000,fp);
      id[i] = i;
      mx[i] = -m200;
      if(isat) mx[i] = -mpeak;
      if(!isat) mhost[i] = m200;
    }

  sort2dbl(nhalo, mx, id);
  for(i=1;i<=nhalo;++i)
    mx[i] = -mx[i];

  deconvolved_smf();

  for(i=1;i<=nbins;++i)
    xbias[i] = hcnt[i] = mbin[i] = 0;
  for(j=1;j<=nmass;++j)
    if(mass[j]>=12.7)break;
  jstart = j;
  ncumu = decon[jstart]*dlogm;
  dlogm = (mass[2] - mass[1]);
  fprintf(stderr,"%d %e %e %e\n",jstart,dlogm,mass[2],mass[1]);
  for(i=1;i<=nhalo;++i)
    {
      density = i/volume;
      while(ncumu<density)
	{
	  jstart--;
	  ncumu += decon[jstart]*dlogm;
	}
      if(jstart<0)break;

      mgal = gasdev(&IDUM)*wpl.a[6]+(mass[jstart]);
      k = id[i];
      printf("BOO %d %d %e %e %e %e %e %e %f %f %f\n",
	     i,jstart,ncumu,density,mx[i],mass[jstart],mgal,decon[jstart],x[k],y[k],z[k]);

      ibin = (int)((mgal-11.25)/0.1) + 1;
      if(ibin<1)continue;
      if(mgal>=11.95)ibin = nbins;
      xbias[ibin] += bias_interp(mhost[k],-1);
      hcnt[ibin]++;
      mbin[ibin]+=mx[i];
    }
  for(i=1;i<=nbins;++i) {
    printf("BIASx %e %e %.0f\n",(xbias[i]/hcnt[i]),mbin[i]/hcnt[i],hcnt[i]);
    shambias[i] = xbias[i]/hcnt[i]; }

  free_dvector(mx,1,nhalo);
  free_dvector(mhost,1,nhalo);
  free_ivector(id,1,nhalo);
  free_dvector(y,1,nhalo);
  free_dvector(x,1,nhalo);
  free_dvector(z,1,nhalo);


}

void shambias_loop()
{
  char fname[1090];
  int i,j;
  FILE *fp;

  wpl.satfac = 1.0;
  for (i=1;i<=30;i++)
    {
      //wpl.a[6] = 0.26; //i/100.0;
      wpl.a[6] = i/100.0;

      //wpl.satfac = i*1.0;

      shmr_minimization();
      bias_from_simulation_shmr();
      sprintf(fname,"bias_sham_%.2f_%.1f.dat",wpl.a[6],wpl.satfac);
      fp = fopen(fname,"w");
      for(j=1;j<=8;++j)
	fprintf(fp,"%e %e\n",shambias[j],shambias2[j]);
      fclose(fp);
    }
  exit(0);
}

void bias_from_simulation_shmr()
{
  static int icall=25;
  int nhalo, i,j,k, isat, nbins=8, ibin;
  double xbias[10], hcnt[10], fac, galcnt=0, satcnt=0;
  float m200, m, mpeak, mhost, mgal, mstar, *mx, mlo, mhi, ngal, cbias, vpeak;
  FILE *fp, *fp1;
  char aa[1000], fname[100];
  unsigned long IDUM = 555;
  float xx[10];

  fp = openfile("../halosubs_multidark_z0.41_D200.dat");
  nhalo = filesize(fp);

  sprintf(fname,"shamfile_%d.dat",icall);
  icall++;
  fp1 = fopen(fname,"w");

  for(i=1;i<=nbins;++i)
    xbias[i] = hcnt[i] = 0;
  for(i=1;i<=nhalo;++i)
    {
      fscanf(fp,"%e %e %e %e %d",&m200, &mhost, &mpeak, &vpeak, &isat);
      for(j=0;j<6;++j)fscanf(fp,"%f",&xx[j]);
      //fgets(aa,1000,fp);
      // get random galaxy for this halo
      m = m200;
      fac = 1;
      if(isat) { m = mpeak; fac = wpl.satfac; }
      mstar = log10(ms_to_mhalo_inversion(m));
      mgal = gasdev(&IDUM)*wpl.a[6]+(mstar);
      fprintf(fp1,"%e %e %e %e %e %e %e %e %d\n",xx[0],xx[1],xx[2],xx[3],xx[4],xx[5],mgal,m,isat);
      ibin = (int)((mgal-11.25)/0.1) + 1;
      if(ibin<1)continue;
      if(mgal>=11.95)ibin = nbins;
      xbias[ibin] += bias_interp(mhost,-1)*fac;
      hcnt[ibin]+=fac;
      galcnt++;
      if(isat) satcnt += fac;
    }
  printf("SATFRAC: %f\n",satcnt*1./galcnt);

  for(i=1;i<=nbins;++i)
    {
      mlo = pow(10.0,11.3+i/10.0-0.1/2); // these are the bin limits (don't use wpl here since Ncen is not called)
      mhi = pow(10.0,11.3+i/10.0+0.1/2);
      if(i==nbins)mhi=pow(10.0,12.9);
      set_up_hod_for_shmr(mlo,mhi,wpl.a);
      ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;      
      shambias[i] = xbias[i]/hcnt[i];
      shambias2[i] = cbias;
      printf("BIASx %e\n",xbias[i]/hcnt[i]);
    }
  fclose(fp);
  fclose(fp1);
}
