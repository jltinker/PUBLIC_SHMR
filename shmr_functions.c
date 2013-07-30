#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* Local functions.
 */
double set_low_mass_shmr(void);
double func_low_mass_shmr(double m);

/* local globals
 */
double ncenhi_g32;
double LO_FRAC = 1.0E-3;

//-----------------------------------------------------------
//  STELLAR MASS TO HALO MASS
//  ms is stellar mass in LINEAR UNITS (non log units)
//-----------------------------------------------------------
//
// input: stellar masses in M_\odot [no h^-2]
// output: halo masses in M_\odot/h, unless otherwise specifed 
double ms_to_mhalo(double ms, double *a)
{
  double res, mhalo_norm, mstellar_norm,beta,delta,gamma,x;

  mhalo_norm     = a[1];          // M1 (in LOG10 units)   
  mstellar_norm  = pow(10,a[2]);  // M*,0 (convert to linear units ...)
  beta           = a[3];          // beta
  delta          = a[4];          // delta
  gamma          = a[5];          // gamma

  x=ms/mstellar_norm;             // in lin units

  // This is in log10
  res = mhalo_norm+(beta*log10(x))+(pow(x,delta)/(1+pow(x,-gamma)))-0.5;

  // Now convert to linear units
  if(HINVERSE_UNITS)
    res = pow(10.0,res)*HUBBLE;
  else
    res = pow(10.0,res);
  return(res);
}



/* This is the red central fraction as a function of halo mass.
 * This is implemented as a spline function interpolated between 6.3e10
 * and 10^14.
 * This code does not assume that the peak galaxy mass at a given halo mass
 * is the same between the reds and blues.
 * ---
 * The values of fredcen are stored in a[16]-a[20] assuming alpha NOT FREE
 * ---
 * Will interpolate in log Mass but in linear of fred. 
 */
// NB!! should i remove the a-vector assumption and simply use the mass vector itself?
// NB!! make the length and the mass limits user-defined!!!!!!!!!!!!
double red_central_fraction(double mass, double *a)
{
  int i,j,k;
  static double *mm, *yy, *aa, aprev=20;
  static int flag=1, n=5, nbuf=22;
  double fred;

  if(COLOR==0)return 0;

  if(flag)
    {
      flag = 0;
      aa = dvector(1,n);
      mm = dvector(1,n);
      for(i=1;i<=n;++i)
	mm[i] = log(1.0e14/6.3e10)/(n-1.0)*(i-1) + log(6.3e10);
      yy = dvector(1,n);
    }
  if(wpl.reset_fred)
    {
      wpl.reset_fred = 0;
      for(i=1;i<=n;++i)
	aa[i] = a[i+nbuf];
      aa[1] = pow(10.0,aa[1]);
      aa[2] = pow(10.0,aa[2]);
      spline(mm,aa,n,1.0E+30,1.0E+30,yy);
    }

  for(i=110;i<=-150;++i)
    {
      splint(mm,aa,yy,n,i/10.0*log(10),&fred);
      fprintf(stdout,"%d %e\n",i,fred);
    }
  for(i=1;i<=-n;++i)
    printf("%e %e %e\n",mm[i],aa[i],yy[i]);

  //if(log(mass)>mm[n])return aa[n];

  splint(mm,aa,yy,n,log(mass),&fred);

  if(fred>1)fred=1;
  if(fred<0)fred=1.0E-10;
  
  return fred;

}

/* qromo called function: Returns the number density of 
 * galaxies within a BIN of stellar mass--> 
 * ie, it will use the cen/sat funcs below
 * NB! This should be made redundant by the fact that Nsat and Ncen 
 * now always call the xigm versions.
 */
double func_galaxy_density_xigm(double m)
{
  m = exp(m);
  //printf(">> %e %e %e\n",m,Ncen_xigm(m),Nsat_xigm(m));
  return (Ncen_xigm(m) + Nsat_xigm(m))*m*dndM_interp(m);
}

/* For a given bin in stellar mass, this function takes the
 * limits of the bin as input, as well as the vector that 
 * corresponds to the SHMR parameter set (11 parameters)
 * 
 * NB! Should I replace BLUE_FLAG and have the code always submit
 * the proper a-vector?
 */
int set_up_hod_for_shmr(double mlo, double mhi, double *a)
{
  int i, ibuf;
  double dmsdmh, m1, m2, mm, mass;

  ibuf = SHMR_PARAMS*(1-BLUE_FLAG);

  wpl.mstar_upper = log10(mhi);
  wpl.mstar_lower = log10(mlo);

  if(CENTRALS_ONLY)
    {
      wpl.mstar = log10(mlo);                 // This is set to retabulate Ncen_CLF for set_low_mass (!only set for lower bin!)
      HOD.M_min = ms_to_mhalo(mlo,&(a[ibuf]));  // convert to h^-1 Msol
      HOD.M_low = set_low_mass_shmr();	  
      return 0;
    }

  // -- construct HOD for low-mass limit --
  wpl.mstar = log10(mlo);                 // This is set to retabulate Ncen_CLF for set_low_mass (!only set for lower bin!)
  HOD.M_min = ms_to_mhalo(mlo,&(a[ibuf]));  // convert to h^-1 Msol
  HOD.MaxCen = red_central_fraction(HOD.M_min,a);
  if(BLUE_FLAG)HOD.MaxCen=1-HOD.MaxCen;

  HOD.M1    = pow((HOD.M_min/1e12),a[10+ibuf])*a[8+ibuf]*1e12;   //linear function for log(M1) and log(Mmin) (note:M1 is Msat)
  if(SHMR_PARAMS>11) {
    HOD.M1 *= 1.0/(1+pow(HOD.M_min/a[12+ibuf],a[13+ibuf]));
    HOD.M1 *= exp(-a[14+ibuf]/HOD.M_min);
  }
  HOD.alpha = a[11+ibuf];
  if(VARIABLE_ALPHA) {
    if(mlo<6.0E10) HOD.alpha = a[11+ibuf]*pow(mlo/6.0E10,a[15+ibuf]);
    else HOD.alpha = a[11+ibuf]*pow(mlo/6.0E10,a[16+ibuf]);
  }
  HOD.M_cut = pow((HOD.M_min/1e12),a[9+ibuf])*a[7+ibuf]*1e12;     //linear function for log(Mcut) and log(Mmin)

  HOD.M_low = set_low_mass_shmr();	  
  if(HOD.M_low<0)return -1;
   
  // -- construct HOD for hi-mass limit --
  // -- use wpl here and not HOD --
  wpl.mmin = ms_to_mhalo(mhi,&(a[ibuf])); // convert to h^-1 Msol
  wpl.maxcen = red_central_fraction(wpl.mmin,a);
  if(BLUE_FLAG)wpl.maxcen=1-wpl.maxcen;

  wpl.m1    = pow((wpl.mmin/1e12),a[10+ibuf])*a[8+ibuf]*1e12;   //linear function for log(M1) and log(Mmin) (note:M1 is Msat)
  if(SHMR_PARAMS>11) {
    wpl.m1 *= 1.0/(1+pow(wpl.mmin/a[12+ibuf],a[13+ibuf]));
    wpl.m1 *= exp(-a[14+ibuf]/wpl.mmin);
  }
  HOD.alpha = a[11+ibuf];
  if(VARIABLE_ALPHA) {
    if(mhi<6.0E10) HOD.alpha = a[11+ibuf]*pow(mhi/6.0E10,a[15+ibuf]);
    else HOD.alpha = a[11+ibuf]*pow(mhi/6.0E10,a[16+ibuf]);
  }
  wpl.mcut = pow((wpl.mmin/1e12),a[9+ibuf])*a[7+ibuf]*1e12;     //linear function for log(Mcut) and log(Mmin)

  wpl.alpha = HOD.alpha;

  //-----

  if(HOD.M_low > HOD.M_min)
    fprintf(stderr,"WARNING: M_low>M_min for Mstar_bin=[%e,%e] %e %e\n",mlo,mhi,HOD.M_low,HOD.M_min);
  
  muh(0);
  //GALAXY_DENSITY = qromo(func_galaxy_density_xigm,log(HOD.M_low),log(HOD.M_max),midpnt); 
  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt); 
  muh(1);
  //GALAXY_BIAS = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  BETA = pow(OMEGA_M,GAMMA)/GALAXY_BIAS;
  muh(2);

  // let's trap the error before then
  if(ERROR_FLAG==1)
    {
      printf("ERROR TRAPPED in setup\n");
      ERROR_FLAG = 0;
      return -1;
    }

  if(ERROR_FLAG == 1)
    {
      printf("\n 1) STOPPED IN set up hod for lensing %e %e\n",wpl.mstar_lower, wpl.mstar_upper); 
      printf("%e %e %e %e %e %e %e %e\n",HOD.M_low,HOD.M_max,HOD.M1, wpl.m1, mlo, mhi, GALAXY_DENSITY, HOD.M_min);
      printf("PARAMS 1 1 ");
      printf("TEST %e\n",qtrap(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),1.0E-5));
      for(i=1;i<=wpl.ncf;++i)
	printf("%.8e ",wpl.a[i]);
      printf("\n");
      for(i=110;i<=150;++i)
	{
	  mm = pow(10.0,i/10.0);
	  printf("HOD %e %e %e\n",mm,Ncen_xigm(mm),Nsat_xigm(mm));
	}
      exit(0);
    }

  return 0;
}


/* Returns N_cen for a stellar mass BIN.
 * NB -- assumes that the HOD currrently describes the LOWER MASS LIMIT of the bin.
 * Mass is h^-1
 */
double Ncen_xigm(double m)
{
  static double mstar_prev=-1;
  static double *mh, *ncen, *zz;
  static int n=0, niter=0;

  int i, ibuf;
  double mlo, mhi, dlogm, a, ms1, ms2, ncen1, ncen2, max;

  ibuf = (1-BLUE_FLAG)*SHMR_PARAMS;

  if(!n)
    {
      n=200;
      mh = dvector(1,n);
      ncen = dvector(1,n);
      zz = dvector(1,n);
    }

  if(mstar_prev != wpl.mstar_lower) // retabulates on lower mass limit
    {
      mstar_prev = wpl.mstar_lower;
      niter++;

      mlo = log(1.0E8);
      mhi = log(1.0E+16);
      dlogm = (mhi-mlo)/(n-1);

      ms1 = wpl.mstar_upper; // assuming that wpl.mstar is in log10 already!!
      ms2 = wpl.mstar_lower; // assuming that wpl.mstar is in log10 already!!

      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-1)*dlogm + mlo);
	  ncen1 = 0.5*(1 - erf((ms1 - log10(ms_to_mhalo_inversion(mh[i])))/(ROOT2*wpl.a[6+ibuf])));
	  ncen2 = 0.5*(1 - erf((ms2 - log10(ms_to_mhalo_inversion(mh[i])))/(ROOT2*wpl.a[6+ibuf])));

	  max = red_central_fraction(mh[i],wpl.a);
	  if(BLUE_FLAG) max = 1-max;
	  ncen[i] = (ncen2 - ncen1)*max;

	  //printf("NCEN %d %e %e %e %e %e %e \n",niter,mh[i],ncen[i],ncen1,ncen2,wpl.a[6+ibuf],log10(ms_to_mhalo_inversion(mh[i])));

	  //Alexie added condition (for Nebula)
	  //Needed to change to -18 here for nebula
	  if(ncen[i] < -1.0E-18){
	    printf("NCEN : %e %e %e %e %e \n",mh[i],ncen1,ncen2,wpl.a[6],log10(ms_to_mhalo_inversion(mh[i])));
	    printf("NCEN : %e %e \n",ms1,ms2);
	    printf("ERROR: negative Ncen in [%f, %f] for M= %e\n",ms2,ms1,m);
	    exit(0);
	  }

	  mh[i] = log(mh[i]);
	  ncen[i] = log(ncen[i] + 1.0E-18);
	}
      spline(mh,ncen,n,1.0E+30,1.0E+30,zz);
    }
  
  splint(mh,ncen,zz,n,log(m),&a);
  return exp(a); 
}

/* This function sets a lookup table for the inversion of the
 * mstar-to-mhalo relation.
 * ----------------------------
 * NB-- for color-defined samples, the MaxCen simply is a multiplicative
 * factor at the end of calculating Ncen(M), so it doesn't need to be 
 * included here.
 */
double ms_to_mhalo_inversion(double mass)
{
  static int n=1000, flag=1;
  static double *ms, *mh, *zz;

  int i, ibuf;
  double mlo, mhi, dlogm, a;


  if(flag || wpl.reset_inversion)
    {
      ibuf = SHMR_PARAMS*(1-BLUE_FLAG);
      //     printf("RESETTING INVERSION: %d\n",ibuf);

      if(flag)
	{
	  ms = dvector(1,n);
	  mh = dvector(1,n);
	  zz = dvector(1,n);
	  flag = 0;
	}
      wpl.reset_inversion = 0;
      
      mlo = log(1.0E-7);
      mhi = log(5.0E12);
      dlogm = (mhi-mlo)/(n-1);

      for(i=1;i<=n;++i)
	{
	  ms[i] = ((i-1)*dlogm + mlo);
	  mh[i] = log(ms_to_mhalo(exp(ms[i]),&(wpl.a[ibuf])));
	  if(i>1)
	    if(mh[i] <= mh[i-1])
	      printf("ERROR: non-monotonic increase in halo mass for inversion!!! %d %e %e\n\n",BLUE_FLAG,mh[i],ms[i-1]);
	}
      spline(mh,ms,n,1.0E+30,1.0E+30,zz);
    }
  splint(mh,ms,zz,n,log(mass),&a);
  return exp(a);

}

/* Returns N_sat for a stellar mass BIN.
 * NB -- assumes that the HOD currrently describes the LOWER MASS LIMIT of the bin.
 * Mass is h^-1.
 * 
 * NB!! -- Note that there is no variation in the functional form of N_sat: 
 * you're stuck with the one written down here, unless you want to alter the code.
 */
double Nsat_xigm(double m)
{

  double mmin, m1, m2, mcut, slogm, mlow, nsat_mlo, nsat_mhi, nsat, dmsdmh, alpha;

  // replacing this with actual function
  nsat_mlo = pow(m/HOD.M1,HOD.alpha)*exp(-(HOD.M_cut+HOD.M_min)/m); // Alexie changed from Mmin to Mcut here

  // put HOD in temp space
  mmin = HOD.M_min;
  m1 = HOD.M1;
  mcut = HOD.M_cut;
  slogm = HOD.sigma_logM;
  mlow = HOD.M_low;
  alpha = HOD.alpha;

  // transfer HOD for upper limit
  HOD.M_min = wpl.mmin;
  HOD.M1 = wpl.m1;
  HOD.M_cut = wpl.mcut;
  HOD.M_low = wpl.mlow;
  HOD.sigma_logM = wpl.slogm;
  HOD.alpha = wpl.alpha;

  // replace this with actual function
  nsat_mhi = pow(m/HOD.M1,HOD.alpha)*exp(-(HOD.M_cut+HOD.M_min)/m); // Alexie changed from Mmin to Mcut here

  // get difference (+ error check)
  nsat = nsat_mlo - nsat_mhi;

  if(isnan(nsat))
    {
      printf("NSAT %e %e %e %e\n",wpl.mstar_lower, wpl.mstar_upper, HOD.M_min, mmin);
      printf("NSAT %e %e %e %e %e %e %e %e\n",m,nsat,HOD.M1,m1,HOD.M_cut,mcut,nsat_mlo, nsat_mhi);
      printf("NSAT %e %e %e %e %e %e %e %e %e %e\n",
	     m,nsat,m1,mcut,alpha,nsat_mlo, HOD.M1,HOD.M_cut,HOD.alpha,nsat_mhi);
      exit(0);
    }
  if(nsat<0)nsat = 0;

  //return the HOD to correct place.
  HOD.M_min = mmin;
  HOD.M1 = m1;
  HOD.M_cut = mcut;
  HOD.M_low = mlow;
  HOD.sigma_logM = slogm;
  HOD.alpha = alpha;

  return nsat;
}


/* This is to find M_low for a given SHMR set and stellar mass threshold.
 * This is a souce of significant consternation in terms of error trapping,
 * so we'll just make this a brute-force calculation.
 */
double set_low_mass_shmr()
{
  double ncenlo, ncenhi, mlo;
  int n=0;

  if(HOD.M_min>HOD.M_max) 
    {
      ncenhi = N_cen(HOD.M_max);
      mlo = HOD.M_max;
    }
  else 
    { 
      ncenhi = N_cen(HOD.M_min);
      mlo = HOD.M_min;
    }
  if(ncenhi<1.0E-17)return -1;
  ncenlo = N_cen(mlo);
  //printf("LOW %e %e %e %e\n",mlo,ncenlo,ncenhi,HOD.M_min);
  while(ncenlo/ncenhi>LO_FRAC) {  
    mlo /= 3.0; 
    ncenlo = N_cen(mlo);     
    //printf("LOW %e %e %e\n",mlo,ncenlo,ncenhi);
    n++;
    //if(n==10)    exit(0);
  }
  ncenhi_g32 = ncenhi;
  //printf("%e %e %e %e\n",HOD.M_min,mlo,N_cen(mlo)/ncenhi,N_cen(mlo*3)/ncenhi);
  return zbrent(func_low_mass_shmr,mlo,mlo*3,1.0E-5);
}
double func_low_mass_shmr(double m)
{
  return N_cen(m)/ncenhi_g32/LO_FRAC - 1;
}
