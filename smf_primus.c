#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"


double boss_bias_function(double logmass, double redshift, double blue_flag);
double boss_redshift_evolution();
void galaxies_from_halos(double mass_halo, double z_real, double rad, int *ngal, double *ra_gal, double *dec_gal, 
			 double *rad_gal, double *vrad_gal);

double chi2_smf_primus()
{
  static int n=0, npr, npb, n_ratio, nr, k, istart_b, istart_r, nb, nbias;
  static double *mstar, *nstar, *estar, **covar, *eval, **evect, *bias_eff;
  static double *mstar_b, *nstar_b, *estar_b, **covar_b;
  static double *mstar_r, *nstar_r, *estar_r, **covar_r;
  static double *mstar_ratio, *nstar_ratio, *estar_ratio, **covar_ratio;
  int i, j, i1, ii, ibuf, nlim = 20, eflag, ngalr;
  FILE *fp;
  static double dlog10m, BOSS_ngaltot, BOSS_ngaltot_err, BOSS_nblue;
  static int niter=0;
  double ng1, ng2, mlo, mhi, ngal, ngalprev, dmsdmh, m1, m2, ns1, ns2, **tmp, **tmp2, mass, nsat, chi1, cbias, 
    x, fsat, fsat_prev, ngal_total, cbias_tot, fac;
  double r, rmin, rmax, dlogr;
  int nwp;
  double boss_bias, boss_bias_err, chi2bias=0, g0, chi2alpha;
  char fname[1000], aa[1000], a1;
  double chi2ngal, xmodel[1000], xmodelr[100],xmodelb[100];
  FILE *outfile, *outfile_cov;
  float x1;

  MCMC_OUTPUT = 1;
  HUBBLE_UNITS = 1;
  MAX_STELLAR_MASS = 1.0E+14;
  wp.pi_max = 80;
  

  g0 = growthfactor(REDSHIFT)/growthfactor(0.55);
  g0 *= sqrt(1.2); // to account for the fact that we're fitting redshift-space bias.


  // put a prior on mhalo star
  if(wpl.a[12]<13)return 1.0E+17;
  

  boss_redshift_evolution();

  /* The HOD set up is different for the stellar mass function because
   * only really need to know the total number of central and satellites
   * not the full Ncen=f(M) and Nsat=f(sat) (doesn't call these functions)
   */
  if(!n)
    {
      //fp = openfile("/Users/tinker/cosmo/CONFORMITY/PRIMUS_SMF/smf_starforming_z0_fit.dat");
      fp = openfile(ARGV[4]);
      nb = filesize(fp);

      // check for header
      if(fgetc(fp)=='#')
	{
	  fgets(aa,1000,fp);
	  fgets(aa,1000,fp);
	  nb = nb - 2;
	}
      else
	rewind(fp);

      nlim = nb;
      mstar_b = dvector(1,nlim);
      nstar_b = dvector(1,nlim);
      estar_b = dvector(1,nlim);

      BOSS_ngaltot = 0;
      for(i=1;i<=nb;++i) {
	fscanf(fp,"%lf %lf %lf",&mstar_b[i],&nstar_b[i],&estar_b[i]);
	mstar_b[i] = fabs(mstar_b[i]) - 10;
	estar_b[i] = nstar_b[i]*0.25 + estar_b[i]*5;
	if(nstar_b[i]<=0) estar_b[i] = 1.0E-7;
	BOSS_ngaltot += nstar_b[i]*0.2;
	fgets(aa,1000,fp); }
      BOSS_nblue = BOSS_ngaltot;
      fclose(fp);
      fprintf(stderr,"Read [%d] lines from starforming SMF \n",nb);

      //fp = openfile("/Users/tinker/cosmo/CONFORMITY/PRIMUS_SMF/smf_quenched_z0_fit.dat");
      fp = openfile(ARGV[5]);
      nr = filesize(fp);

      // check for header
      if(fgetc(fp)=='#')
	{
	  fgets(aa,1000,fp);
	  fgets(aa,1000,fp);
	  nr = nr - 2;
	}
      else
	rewind(fp);

      nlim = nr;
      mstar_r = dvector(1,nlim);
      nstar_r = dvector(1,nlim);
      estar_r = dvector(1,nlim);
      
      for(i=1;i<=nr;++i) {
	fscanf(fp,"%lf %lf %lf",&mstar_r[i],&nstar_r[i],&estar_r[i]);
	mstar_r[i] = fabs(mstar_r[i]) - 10;
	estar_r[i] = nstar_r[i]*0.1 + estar_r[i]*5;
	if(nstar_r[i]<=0) estar_r[i] = 1.0E-7;
	BOSS_ngaltot += nstar_r[i]*0.2;
	fgets(aa,1000,fp); }
      fclose(fp);
      fprintf(stderr,"Read [%d] lines from quenched SMF \n",nr);
	  
      // using the tabulated data, determine the binsize
      dlog10m = mstar_b[2] - mstar_b[1];
      fprintf(stderr,"mcmc> SMF binsize: %f\n",dlog10m);
      n = 1;

      // let get the overall total number of galaxies
      BOSS_ngaltot_err = 0.01*BOSS_ngaltot;

      fp = openfile("bias_model_total.dat");
      nbias = filesize(fp);
      bias_eff = dvector(0,nbias-1);
      for(i=0;i<nbias;++i) {
	fscanf(fp,"%f %lf",&x1,&bias_eff[i]); fgets(aa,1000,fp); }
      fclose(fp);
      printf("EFFECIVE BIAS for IZED=%d b=%f\n",IZED,bias_eff[IZED]);

    }
  wpl.magnitude_cutoff = wpl.a[28];
  wpl.cutoff_width = wpl.a[29];
  

  chi2ngal = 0;
  ngal = 0;

  /*** ---- BLUE STELLAR MASS FUNCTION ----- ***/

  BLUE_FLAG = 1;
  wpl.reset_inversion = 1;
  i=set_up_hod_for_shmr(1.0E+10+drand48()*1.0E+10,6.0E10,wpl.a);
  //if we have an issue, return big chi2 value
  if(i==-1){ niter++; return 1.0E+17; }
  
  for(i=1;i<=n;++i)
    xmodel[i] = 0;

  // only do the overall blue gal bias
  chi2bias = 0;
  mlo = pow(10.0,mstar_b[1]-dlog10m/2);
  mhi = pow(10.0,mstar_b[nb]+dlog10m/2);
  if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
  if(HOD.M_low<0)
    chi2ngal = 1.0E+17;

  // only do centrals for the blues
  ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;      


  // only do the overall number density of blue galaxies
  // don't worry about completeness cutoff.
  // use same error as the bias
  for(ng1=0, i=1;i<=nb;++i)
    ng1 += nstar_b[i]*dlog10m;
  x = 0.01*ng1 + 0.01*4.0E-4;
  x = 100000;

  ngal_total = 0;
  cbias_tot = 0;
  for(i=1;i<=nb;++i)
    {

      mlo = pow(10.0,mstar_b[i]-dlog10m/2);
      mhi = pow(10.0,mstar_b[i]+dlog10m/2);

      eflag = ngal  = 0;
      if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
      if(HOD.M_low<0)
	{
	  chi2ngal += 1.0E+17;
	  continue;
	}
      //ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);

      ngalprev = ngal;
      ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      //cbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;      
      ngal = ngal*pow(HUBBLE_UNITS,3.0)/dlog10m; // Galaxy density, convert to h72 units
      ngal = ngal*0.5*(1+erf((mstar_b[i]-wpl.magnitude_cutoff)/wpl.cutoff_width));
      ngal_total += ngal*dlog10m;
      
      // only do centrals for the blues
      ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal*g0;      
      cbias_tot += cbias*ngal*0.5*(1+erf((mstar_b[i]-wpl.magnitude_cutoff)/wpl.cutoff_width));;
      // printf("BIAS %e %e %e %e\n",cbias_tot,cbias,ngal,ngal_total);
      continue;

      // check for monotonically decreasing mass function
      /*
      if(i>n)
	{
	  if(ngalprev<ngal)chi2ngal += 1.0E6; // previous density should be higher
	  break;
	}
      */  
      if(HOD.M_low>HOD.M_max)
	{
	  if(ERROR_FLAG==1)ERROR_FLAG = 0;
	  if(isnan(ngal))ngal = 0;
	}	  


      // convert to log10 (ONLY FOR PRIMUS)
      //ngal = log10(ngal);
      chi2ngal += (ngal-nstar_b[i])*(ngal-nstar_b[i])/(estar_b[i]*estar_b[i]);
      xmodel[i] = ngal;

      if(MCMC_OUTPUT)
	printf("SMFB%d %e %e %e %e %e %f %f\n",niter,mstar_b[i],nstar_b[i],estar_b[i],ngal,
	       ms_to_mhalo(pow(10.0,mstar_b[i]),wpl.a),cbias,boss_bias);

      if(ERROR_FLAG == 1){
	printf("\n STOPPED IN chi2_stellar_mass_function %d\n",i);
	printf("SMF%d %e %e %e %e %e\n",niter,mstar[i],nstar[i],estar[i],ngal,(ns1-ns2)*pow(HUBBLE_UNITS,3.0)/dlog10m);
	printf("%e %e %e %e\n",HOD.M_low, HOD.M_min, HOD.M_max, N_sat(HOD.M_max));
	ERROR_FLAG=0;
	//exit(0);
      }
      //if(ngal<-15)break;
      
    }
  cbias_tot /= ngal_total;
  boss_bias = boss_bias_function(mstar_b[1]+0.5*dlog10m,REDSHIFT,BLUE_FLAG);
  // have the error on blue bias depend on the bias fraction
  boss_bias_err = 0.03*boss_bias;
  if(BOSS_nblue/BOSS_ngaltot<0.1) boss_bias_err = 0.1*boss_bias;
  if(BOSS_nblue/BOSS_ngaltot<0.03) boss_bias_err = 0.5*boss_bias;
  
  chi2bias += (boss_bias-cbias_tot)*(boss_bias-cbias_tot)/(boss_bias_err*boss_bias_err);
  printf("CHIBBIAS%d %e %f %f %f %e %e %e\n",niter,chi2bias,boss_bias, cbias_tot, boss_bias_err, BOSS_nblue, ngal_total, BOSS_ngaltot);
  cbias_tot *= ngal_total;


  if(MCMC_OUTPUT)
    printf("SMFB%d %e %e %e %e %e %f %f\n",niter,mstar_b[1],ng1,x,ngal_total,
	   ms_to_mhalo(pow(10.0,mstar_b[1]),wpl.a),cbias_tot,boss_bias);
  // don't do the chi2ngal until the total
  // unless the blue fraction is above 15%
  if(BOSS_nblue/BOSS_ngaltot>0.15)
    chi2ngal += (ngal_total-BOSS_nblue)*(ngal_total-BOSS_nblue)/(BOSS_nblue*BOSS_nblue*0.07*0.07);


  /*** ---- RED STELLAR MASS FUNCTION ----- ***/
 
  BLUE_FLAG = 0;
  ibuf = SHMR_PARAMS*(1-BLUE_FLAG);
  wpl.reset_inversion = 1;
  i=set_up_hod_for_shmr(1.0E+10+drand48()*1.0E+10,3.0E10,wpl.a);
  //if we have an issue, return big chi2 value
  if(i==-1){ niter++; return 1.0E+17; }

  for(i=1;i<=n;++i)
    xmodel[i] = 0;


  //ngal_total = 0;
  fsat = 1;
  for(i=1;i<=nr;++i)
    {

      mlo = pow(10.0,mstar_r[i]-dlog10m/2);
      mhi = pow(10.0,mstar_r[i]+dlog10m/2);
      
      eflag = ngal = cbias = 0;
      if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
      //ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      if(HOD.M_low<0)
	{
	  printf("%e %e %e %e\n",log10(mlo)+10,log10(mhi)+10,HOD.M_min,nstar_r[i]);
	  chi2ngal += 1.0E+17;
	  continue;
	}      

      ngalprev = ngal;
      ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      cbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal*g0;

      //ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      //cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;      
      ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      
      fsat_prev = fsat;
      fsat = 0;
      if(fsat_prev>0.01)
	fsat = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;
      
      // enforce monotonicity of fsat (within 1%)
      if(fsat>fsat_prev*1.1 && nstar_r[i]>0)
	chi2ngal+=1.0E+17;
      // enforce upper limit on fsat
      if(fsat>0.4 && nstar_r[i]>0)
	chi2ngal+=1.0E+17;

      

      ngal = ngal*pow(HUBBLE_UNITS,3.0)/dlog10m; // Galaxy density, convert to h72 units
      fac = 0.5*(1+erf((mstar_r[i]-wpl.magnitude_cutoff)/(wpl.cutoff_width/1)));//  + 0.04;
      //fac = fac*0.99 + 0.01;
      ngal = ngal*fac;
      ng1 = ng2 = ns1 = ns2 = 0;
      ngal_total += ngal*dlog10m;
      cbias_tot += cbias*ngal*dlog10m;

      // check for monotonically decreasing mass function
      /*
      if(i>n)
	{
	  if(ngalprev<ngal)chi2ngal += 1.0E6; // previous density should be higher
	  break;
	}
      */  
      if(HOD.M_low>HOD.M_max)
	{
	  if(ERROR_FLAG==1)ERROR_FLAG = 0;
	  if(isnan(ngal))ngal = 0;
	}	  

      // convert to log10 (ONLY FOR PRIMUS)
      //ngal = log10(ngal);
      chi2ngal += (ngal-nstar_r[i])*(ngal-nstar_r[i])/(estar_r[i]*estar_r[i]);
      xmodel[i] = ngal;
      
      boss_bias = boss_bias_function(mstar_r[i]+0.5*dlog10m,REDSHIFT,BLUE_FLAG);            
      boss_bias_err = 0.02*boss_bias;      
      if(i==1)
	boss_bias_err = 0.01*boss_bias;
      // don't worry about bias>3
      if(boss_bias<=3.3 && nstar_r[i]>0)
	chi2bias += (boss_bias-cbias)*(boss_bias-cbias)/(boss_bias_err*boss_bias_err);

      if(MCMC_OUTPUT)
	printf("SMFR%d %e %e %e %e %e %f %f %f\n",niter,mstar_r[i],nstar_r[i],estar_r[i],
	       ngal,ms_to_mhalo(pow(10.0,mstar_r[i]),&(wpl.a[ibuf])),cbias,boss_bias,fsat);

      // lets output the HOD
      /*
      sprintf(fname,"hod_%d.dat",i);
      fp = fopen(fname,"w");
      rmin = 12.0;
      rmax = 16.0;
      nwp = 100;
      dlogr = (rmax-rmin)/(nwp-1);
      for(j=1;j<=nwp;++j)
	{
	  r = pow(10.0,j*dlogr+rmin);
	  fprintf(fp,"%e %e %e\n",r,N_cen(r), N_sat(r));
	}
      fclose(fp);
      */      

      // let's output the correlation functions	
      /*
      sprintf(fname,"wp_%d.dat",i);
      fp = fopen(fname,"w");
      rmin = 0.1;
      rmax = 30;
      nwp = 50;
      dlogr = log(rmax/rmin)/(nwp-1);
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      for(j=1;j<=nwp;++j)
	{
	  r = exp((j-1)*dlogr)*rmin;
	  fprintf(fp,"%e %e\n",r,projected_xi(r));
	}
      fclose(fp);
      */

      if(ERROR_FLAG == 1){
	printf("\n STOPPED IN chi2_stellar_mass_function %d\n",i);
	printf("SMF%d %e %e %e %e %e\n",niter,mstar_r[i],nstar_r[i],estar_r[i],ngal,(ns1-ns2)*pow(HUBBLE_UNITS,3.0)/dlog10m);
	printf("%e %e %e %e\n",HOD.M_low, HOD.M_min, HOD.M_max, N_sat(HOD.M_max));
	ERROR_FLAG=0;
      }
      // if(ngalr<-15)break;
      
    }

  // do an overall check on total red ngal
  /*
  for(ng1=0, i=1;i<=nb;++i)
    ng1 += nstar_r[i]*dlog10m;
  x = 0.01*ng1;
  chi2ngal += (ngal_total-ng1)*(ngal_total-ng1)/(x*x);
  */
  // have chi2 be overall total
  chi2ngal += (ngal_total-BOSS_ngaltot)*(ngal_total-BOSS_ngaltot)/(BOSS_ngaltot_err*BOSS_ngaltot_err);
  
  // match the overall bias as well
  cbias_tot /= ngal_total;
  chi2bias += (cbias_tot-bias_eff[IZED])*(cbias_tot-bias_eff[IZED])/(0.015*0.015*bias_eff[IZED]*bias_eff[IZED]);
  printf("BIASTOT%d %e %e %e\n",niter,bias_eff[IZED],cbias_tot,
	 (cbias_tot-bias_eff[IZED])*(cbias_tot-bias_eff[IZED])/(0.015*0.015*bias_eff[IZED]*bias_eff[IZED]));

  /* let's check on monotonicty of red_cen_frac
   */
  if(wpl.a[24]<wpl.a[23])chi2ngal+=1.0E+17;
  if(log10(wpl.a[25])<wpl.a[24])chi2ngal+=1.0E+17;
  if(wpl.a[26]<wpl.a[25])chi2ngal+=1.0E+17;
  if(wpl.a[27]<wpl.a[26])chi2ngal+=1.0E+17;

  /****---- do the red/blue satellite fraction -----***/

  // put the prior on asat
  chi2alpha = 0;
  if(HOD.alpha>1.5) { chi2alpha = (HOD.alpha-1.5)*(HOD.alpha-1.5)/0.1/0.1; }

  if(MCMC_OUTPUT)
    printf("SMFCHI%d %e %e %e %e\n",niter,chi2ngal,chi2bias,chi2alpha,chi2bias+chi2ngal+chi2alpha);
  niter++;
  //exit(0);
  return chi2ngal+chi2bias+chi2alpha;
}

double boss_bias_function(double logmass, double redshift, double blue_flag)
{
  double bb, ll;
  static int prev_cosmo=-100;
  static double g0;

  if(blue_flag) return 1.86;

  ll = pow(10.0,0.4*logmass)/pow(10.0,0.4*11.0);
  bb  = 2.0*(1.04548 + 0.004197*pow(ll,3.14446));

  return bb;
}

double boss_redshift_evolution()
{
  FILE *fp, *fpt, *fph, *fphb;
  int i, j, i1, j1, nmag, eflag, nhalo, k;
  double redshift, maghi, maglo, dmag, ngaltot, biastot, mlo, mhi, mstar, ngal, cbias, fac, g0, bred;
  double a[100];
  char fname[1000];

  double ra_gal, dec_gal, ncen, nsat;
  
  double masslo, masshi, dlogm, mass, ncentot, nsattot;

  masslo = 1.0E12;
  masshi = 5.0E+15;
  nhalo = 1000;
  dlogm = log(masshi/masslo)/(nhalo-1);

  maglo = 10.1;
  maghi = 13.9;
  // nmag = (maghi-maglo)/dmag;
  nmag =38;
  dmag = (maghi-maglo)/nmag;
  //printf("DMAG %f %f %f %d\n",dmag, maghi, maglo, nmag);
  
  REDSHIFT = 0.55;
  RESET_COSMOLOGY++;
  SIGMA_8 = 0.8*growthfactor(REDSHIFT);

  fpt = fopen("CLF_lookup_table.dat","w");
  fph = fopen("HOD_lookup_table_red.dat","w");
  fphb = fopen("HOD_lookup_table_blue.dat","w");

  fprintf(fph,"%d %f %f\n",28,0.21,0.73);
  fprintf(fph,"%d %e %e\n",nhalo,masslo,masshi);
  fprintf(fphb,"%d %f %f\n",28,0.21,0.73);
  fprintf(fphb,"%d %e %e\n",nhalo,masslo,masshi);
  fprintf(fpt,"%d %f %f\n",28,.21,.73);
  fprintf(fpt,"%d %e %e\n",nhalo,masslo,masshi);
  fprintf(fpt,"%d %f %f\n",nmag,maglo+dmag/2,maghi-dmag/2);

  for(j=21;j<=21;++j)
    {
      sprintf(fname,"best.%d",j);
      fp = openfile(fname);
      fscanf(fp,"%d %d",&i1,&j1);
      for(i=1;i<=wpl.ncf;++i)
	fscanf(fp,"%lf",&a[i]);
      fclose(fp);
      for(i=1;i<=wpl.ncf;++i)wpl.a[i] = a[i];

      REDSHIFT = j*0.02 + 0.21;
      SIGMA_8 = 0.8*growthfactor(REDSHIFT);
      RESET_COSMOLOGY++;
      
      g0 = growthfactor(REDSHIFT)/growthfactor(0.55)*sqrt(1.2);


      wpl.magnitude_cutoff = wpl.a[28];
      wpl.cutoff_width = wpl.a[29];
      wpl.reset_fred = 1;

      BLUE_FLAG = 0;
      wpl.reset_inversion = 1;

      // lets create the lookup table of CLFs as function of halo mass, redshift
      for(i=1;i<=nhalo;++i)
	{
	  mass = exp((i-1)*dlogm)*masslo;
	  ncentot = nsattot = 0;
	  for(k=1;k<=nmag;++k)
	    {
	      mstar = (k-0.5)*dmag + maglo;
	      mlo = pow(10.0,mstar-dmag/2);
	      mhi = pow(10.0,mstar+dmag/2);
	      fac = 0.5*(1+erf((mstar-wpl.magnitude_cutoff)/(wpl.cutoff_width/1)));
	      if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
	      fprintf(fpt,"%.2f %e %e %e %e\n",REDSHIFT, mass, mstar+10, N_cen(mass)*fac,N_sat(mass)*fac);
	      if(mass>HOD.M_low)
		{
		  ncentot += N_cen(mass)*fac;
		  nsattot += N_sat(mass)*fac;
		}
	    }
	  fprintf(fph,"%d %.2f %e %e %e\n",j, REDSHIFT,mass,ncentot,nsattot);
	  //fflush(fpt);
	  //fprintf(stdout,"HOD%d %.2f %e %e %e %e %e\n",j, REDSHIFT,mass,ncentot,nsattot, ncen, nsat);
	}


      // now integrate over the luminosity function      
      ngaltot = 0;
      biastot = 0;
      // do the reds first
      BLUE_FLAG = 0;
      wpl.reset_inversion = 1;
      for(i=1;i<=nmag;++i)
	{	  
	  mstar = (i-0.5)*dmag + maglo;
	  mlo = pow(10.0,mstar-dmag/2);
	  mhi = pow(10.0,mstar+dmag/2);
	  if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
	  fac = 0.5*(1+erf((mstar-wpl.magnitude_cutoff)/(wpl.cutoff_width/1)));
	  ngal = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  cbias = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;
	  if(cbias>5)cbias=5;
	  ngaltot += ngal*fac;
	  biastot += cbias*ngal*fac;
	  printf("SMFR%d %e %e %e %e %e\n",j,mstar,ngal*fac/dmag,cbias*g0,biastot/ngaltot*g0,ms_to_mhalo(pow(10.0,mstar),a));
	  //if(ngal<1.0E-7 && mlo>13)break;
	  
	}
      bred = biastot/ngaltot*g0;


      // do the blue galaxies
      BLUE_FLAG = 1;
      wpl.reset_inversion = 1;

      // lets create the lookup table of CLFs as function of halo mass, redshift
      for(i=1;i<=nhalo;++i)
	{
	  mass = exp((i-1)*dlogm)*masslo;
	  ncentot = nsattot = 0;
	  for(k=1;k<=nmag;++k)
	    {
	      mstar = (k-0.5)*dmag + maglo;
	      mlo = pow(10.0,mstar-dmag/2);
	      mhi = pow(10.0,mstar+dmag/2);
	      fac = 0.5*(1+erf((mstar-wpl.magnitude_cutoff)/(wpl.cutoff_width/1)));
	      if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
	      //fprintf(fpt,"%.2f %e %e %e %e\n",REDSHIFT, mass, mstar+10, N_cen(mass)*fac,N_sat(mass)*fac);
	      if(mass>HOD.M_low)
		{
		  ncentot += N_cen(mass)*fac;
		  nsattot += N_sat(mass)*fac;
		}
	    }
	  fprintf(fphb,"%d %.2f %e %e %e\n",j, REDSHIFT,mass,ncentot,nsattot);
	  //fprintf(stdout,"HOD%d %.2f %e %e %e %e %e\n",j, REDSHIFT,mass,ncentot,nsattot, ncen, nsat);
	}




      for(i=1;i<=nmag;++i)
	{
	  mstar = (i-0.5)*dmag + maglo;
	  mlo = pow(10.0,mstar-dmag/2);
	  mhi = pow(10.0,mstar+dmag/2);
	  fac = 0.5*(1+erf((mstar-wpl.magnitude_cutoff)/(wpl.cutoff_width/1)));
	  if(set_up_hod_for_shmr(mlo,mhi,wpl.a)==-1) eflag = 1;
	  ngal = qromo(func_central_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  cbias = qromo(func_central_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/ngal;
	  if(cbias>5)cbias=5;
	  ngaltot += ngal*fac;
	  biastot += cbias*ngal*fac;
	  printf("SMFB%d %e %e %e %e\n",j,mstar,ngal*fac/dmag,cbias,biastot/ngaltot*g0);
	  //if(ngal<1.0E-7 && mlo>13)break;
	}

      biastot = (biastot/ngaltot);
      printf("BOO %d %e %e %e %e\n",j,(j+0.5)*0.02 + 0.2,ngaltot,biastot*g0,bred);

    }
  exit(0);

  return 0;

}


void galaxies_from_halos(double mass_halo, double z_real, double rad, int *ngal, double *ra_gal, double *dec_gal, 
			 double *rad_gal, double *vrad_gal)
{
  static int nhalo=300, nmag=19, flag=1, nz=26;
  static float **hod_cen, **hod_sat, ***clf_cen, ***clf_sat, dlogm, dz, zlo, zhi, masslo, masshi, *redshift, *hmass;
  FILE *fp;
  int i,j,i1, iz, im;
  float ix, nc1,nc2,nc,ns1,ns2,ns;

  if(flag)
    {
      flag = 0;
      fprintf(stderr,"g_from_h> initializing\n");
      fp = fopen("HOD_lookup_table.dat","r");
      // read in the header information
      fscanf(fp,"%d %f %f",&nz,&zlo,&zhi);
      fscanf(fp,"%d %f %f",&nhalo,&masslo, &masshi);
      dlogm = log(masshi/masslo)/(nhalo-1);
      fprintf(stderr,"g_from_h> initializing: %d %d\n",nz,nhalo);

      hod_cen = matrix(1,nz,1,nhalo);
      hod_sat = matrix(1,nz,1,nhalo);
      redshift = vector(1,nz);
      hmass = vector(1,nhalo);

      fprintf(stderr,"g_from_h> starting read in\n");
      for(i=1;i<=nz;++i)
	{	  
	  for(j=1;j<=nhalo;++j)
	    fscanf(fp,"%d %f %f %f %f",&i1,&redshift[i],&hmass[j],&hod_cen[i][j],&hod_sat[i][j]);
	}
      fprintf(stderr,"g_from_h> done with init\n");
      dz = redshift[2] - redshift[1];
      fclose(fp);
    }
  
  // find the redshift bins
  iz = (z_real - zlo - dz/2)/dz;
  if(iz<=0)iz=1;
  if(iz>=nz)iz=nz-1;

  // find nearest mass bin
  ix = (log(mass_halo)-log(masslo))/dlogm + 1;
  im = ix;
  if(ix-(int)ix > 0.5)im = (int)ix + 1;
  if(im<1)im=1;
  if(im>nhalo)im=nhalo;

  // get the mean cen and interpolate
  nc1 = hod_cen[iz][im];
  nc2 = hod_cen[iz+1][im];
  nc = (nc2-nc1)/dz*(z_real-redshift[iz]) + nc1;
  
  ns1 = hod_sat[iz][im];
  ns2 = hod_sat[iz+1][im];
  ns = (ns2-ns1)/dz*(z_real-redshift[iz]) + ns1;
  
  //print checks
  //printf("> %f %d %f %f\n",z_real,iz,redshift[iz],redshift[iz+1]);
  //printf("> %e %d %f %e %e\n",mass_halo,im,ix,hmass[im],hmass[im+1]);
  //printf("> %e %e %e\n",nc1,nc2,nc);
  //printf("> %e %e %e\n",ns1,ns2,ns);

  *ngal = poisson_deviate(ns);
  if(drand48()>nc)(*ngal)++;
}

