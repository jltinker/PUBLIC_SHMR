#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

/* Internal functions */
double tabulate_xi(double r);                      // Tabulate xi (for the first part of the integral)
double tabulate_xi_2h(double r);                   // Tabulate xi for the 2h term (to speed up the calculation)
void calc_xi(double *r, double *xi, int n);        // Returns xi in linear units and r in linear Mpc
double xi_1h(double r);                            // Calculate xi for one halo
double xi_1h_simple(double r);                     // Testing xi_1h for NFW
double funct_1h(double m);                         // Function that defines 1h term to be integrated
double delta_sigma(double r);                      // Calculates Delta Sigma
double sigma(double r);                            // Calculates Sigma
double funct1(double r);                           // Funct for DS (calculate mean within sphere of radius r)
double xi2sigmabar(double r);                      // From xigm to Sigma_bar
double sigmabar_funct(double r);                   // Function to be integrated to go from xi to Sigma_bar
void two_halo_ds(double r[8], double *ds);         // Calcuates two halo Delta Sigma
void one_halo_ds(double r[8], double *ds);         // Calculates one halo Delta Sigma
void full_mod_ds(double r[8], double *ds);         // Integrates over 1h and 2h together
void ps_ds(double r[8], double ms, double *ds);    // Point source term
double ms_to_mhalo(double ms, double *a);          // Link stellar mass to halo mass
void make_pretty_plot(double ms, int plotnum);     // Output for a nice plot
void sigma_fabian(double ms);         // Sigma for Fabian
void make_pretty_plot_foredo(double ms);
void ds_mod(double r[8], double ms, double *ds, double *res_point_source,  double *res_2h, double *res_1h, double *res_full, double *a);  // Full DS model
void write_text_halo_concentration(void);          // Write text file for get_conc.pro
void write_sigma_pdf(double r[8]);                 // * NEED TO FINISH THIS *


double chi2_lensing_red(double *a);
double chi2_lensing_blue(double *a);

/* test functions */
void test_xi(double n);
void test_delta_sigma(void);

/* Local global variable*/
double alex_trans_R, alex_mmin, alex_mmax, alex_ds_r, alex_sat_r;
int ALEX_RESET_FLAG_XI, ALEX_FLAG_2H_TABULATE, ALEX_XI_OPTION, ALEX_PRETTY_PLOT_FLAG, ALEX_2H_BIN_FLAG, ALEX_USE_2H_TABULATE,ALEX_XI_TEST_OPTION;
int N_LENSING, N_LENSING_DATA_BINS;               // Number of lensing data points

/* Local globale for using covariance matrix
 */
int LENSING_COVAR = 0;

// ALEX_RESET_FLAG_XI    : tabulate x1 again
// ALEX_FLAG_2H_TABULATE : tabulate the xi 2h again
// ALEX_XI_OPTION        : 1= 1 halo, 2= 2 halo, 3= integrate together
// ALEX_PRETTY_PLOT_FLAG : used to make data points for pretty plot (seperate cent and sat terms)
// ALEX_2H_BIN_FLAG      : which sm bin to use for 2h tabulate
// ALEX_USE_2H_TABULATE  : implement the 2h tabulate version (set to 1 at top of program)
// ALEX_XI_TEST_OPTION   : test the calculation of the NFW profile, use with test_chi2_lensing.pro

// note: Jeremy Flag: RESET_FLAG_2H: to reset the calculation of the 2h term

/* Note LENSING_OUTPUT_FLAG is set at the command line */

// AUG 20: changed all the dump = xi_2h_gm(1, alex_mmin, alex_mmax, wpl.a) to set_lensing_hod

// NOTE : RESET_FLAG_2H is made redundant f=by set_up_hod_for_lensing ... ## actually maybe used in two_halo_rspace

// Timing of code :
// All the time bottleneck is in the term xi_2h_gm -> this is 3s if we use the full exclusion model!
// Can also choose to tabulate the 2h term since is does not affect the g-g lensing very much
// After xi_2h_gm: integrating full_ds_model takes 0.8s * 7 bins = 5.6s (actually I think this timing estimate is wrong? need to redo)


/**---------------------------------------------------------
 **
 ** do the chi^2 lensing for the BLUE galaxies
 **
 **--------------------------------------------------------*/

// JPL: replace /Users/alexie/Work/
// /home/asleauth/

double chi2_lensing(double *a)
{
  return chi2_lensing_blue(a) + chi2_lensing_red(a);
}

double chi2_lensing_blue(double *a)
{
  int i,j,test,n,i1;
  char fname[1000], aa[1000];                                 // To read in data
  FILE *fp;                                                   // To read in data
  double z,dump;                                              // Redshift and dump vector
  double **tmp,**tmp2;                                        // Dump matrix for inverting covar matrix
  double t1,t2,t3,t4;                                         // Time code
  double ds_sat_temp,sig_sat_temp,sig_cen_temp,sig_all_temp;  // Calcualte 1h_sat for plotting purposes
  double chi0,chi1,chi2,chi3,chi4,chi5,chi6,chi_tot;          // Chi2 for various stellar mass bins
  static int FIRST_LENSING_RUN=1;                             // Set to 1 and put to 0 after data has been read in (!! set to one here to only read in once)
  static float Z_CUT1,Z_CUT2,Z_CUT3,Z_CUT4;                   // Redshift bins
  static double *xx;                                          // dump vector
  static double sm_lim_0,sm_lim_1,sm_lim_2,sm_lim_3,sm_lim_4,sm_lim_5,sm_lim_6,sm_lim_7;  // Stellar mass bins
  static double sm0,sm1,sm2,sm3,sm4,sm5,sm6;                                              // Mean stellar mass in each bin
  static double sm0_log,sm1_log,sm2_log,sm3_log,sm4_log,sm5_log,sm6_log;                  // Log stellar mass
  static double *r0,*ds0,*err_ds0, *res_ds0, *res_ps0, *res_2h0, *res_1h0, *res_full0;    // Data and fits to data
  static double *r1,*ds1,*err_ds1, *res_ds1, *res_ps1, *res_2h1, *res_1h1, *res_full1;
  static double *r2,*ds2,*err_ds2, *res_ds2, *res_ps2, *res_2h2, *res_1h2, *res_full2;
  static double *r3,*ds3,*err_ds3, *res_ds3, *res_ps3, *res_2h3, *res_1h3, *res_full3;
  static double *r4,*ds4,*err_ds4, *res_ds4, *res_ps4, *res_2h4, *res_1h4, *res_full4;
  static double *r5,*ds5,*err_ds5, *res_ds5, *res_ps5, *res_2h5, *res_1h5, *res_full5;
  static double *r6,*ds6,*err_ds6, *res_ds6, *res_ps6, *res_2h6, *res_1h6, *res_full6;
  static double **covar0,**covar1,**covar2,**covar3,**covar4,**covar5,**covar6;           // Covar matrix

  if(DONT_FIT_LENSING)
    return(0.0);

  wpl.reset_inversion = 1;
  BLUE_FLAG = 1;
  EXCLUSION=5;                 // HOD exclusion: this calculates 2h for galaxy-matter cross-correlation 
  test=0;                      // Print some test statements
  N_LENSING =8;                // Number of lensing data points
  ALEX_USE_2H_TABULATE=1;      // Implement the tabulated version of 2h
  ALEX_XI_TEST_OPTION=0;       // Use to test delta_sigma calculation
  
  ALEX_PRETTY_PLOT_FLAG =0;    // (*Don't set here*) Set this later to call 1h_cent and 1h_sat

  //printf(" -- NOW ENTERING CHI2_LENSING\n");
  //t1=second();

  // Use this to make a M-c relation for IDL
  //write_text_halo_concentration();

  // Test Halo Mass Function
  /*  fp=fopen("/Users/alexie/Work/HOD/Pro/dndm_code_0.88.txt","w");
  for(i=0;i<=100;i++){
    fprintf(fp,"%lf %lf \n",10+(5.*i/100.0),log10(dndM_interp(pow(10.0,10+(5.*i/100.0)))));
  }
  fclose(fp);
  exit(0);*/

  if(ERROR_FLAG == 1){
    printf("\n STOPPED at the top of chi2_lensing. There is a previously undetected ERROR_FLAG\n");
    exit(0);
  }

  Z_CUT1 = 0.22;               // Redshift bins
  Z_CUT2 = 0.48; 
  Z_CUT3 = 0.74;
  Z_CUT4 = 1.0;
  
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 7;
  }
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 7;
  }
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 6;
  }

  // These are the parameters for the final plot
  // This is set in MCMC LENSING
  if(LENSING_OUTPUT_FLAG){
    // Change N_lensing for finner binning (see the files: finebin.txt)
    N_LENSING = 61; 
  }

  chi_tot = 0;
  chi0 = 0;
  chi1 = 0;
  chi2 = 0;
  chi3 = 0;
  chi4 = 0;
  chi5 = 0;
  chi6 = 0;
   
  // READ IN DATA
  // /home/asleauth/Weak_lensing/Results_gg/SM_PAPER
  // NOTE: these distances are: h=0.7, kpc, physical

  if(FIRST_LENSING_RUN == 1){

    // Stellar mass bins
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      //cut1 = [8.70, 9.20, 9.82, 10.30, 10.64, 10.89, 11.12, 12.0]    
      sm_lim_0 = pow(10,  12.0); // hi
      sm_lim_1 = pow(10,  11.12);
      sm_lim_2 = pow(10,  10.89);
      sm_lim_3 = pow(10,  10.64);
      sm_lim_4 = pow(10,  10.30);
      sm_lim_5 = pow(10,  9.82);
      sm_lim_6 = pow(10,  9.20);
      sm_lim_7 = pow(10,  8.70); //low
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      //cut2 = [9.30, 9.80, 10.30, 10.65, 10.88, 11.05, 11.29, 12.0]
      sm_lim_0 = pow(10,  12.0);
      sm_lim_1 = pow(10,  11.29);
      sm_lim_2 = pow(10,  11.05);
      sm_lim_3 = pow(10,  10.88);
      sm_lim_4 = pow(10,  10.65);
      sm_lim_5 = pow(10,  10.30);
      sm_lim_6 = pow(10,  9.80);
      sm_lim_7 = pow(10,  9.30);
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      // NOTE: only  6 bins here
      //cut3 = [9.80, 10.39, 10.74,10.97,11.16,11.35,12.0]
      sm_lim_0 = pow(10,  12.0);
      sm_lim_1 = pow(10,  11.35);
      sm_lim_2 = pow(10,  11.16);
      sm_lim_3 = pow(10,  10.97);
      sm_lim_4 = pow(10,  10.74);
      sm_lim_5 = pow(10,  10.39);
      sm_lim_6 = pow(10,  9.80);
    }

    xx = dvector(0,N_LENSING-1); //dump vector

    // The various lensing stellar mass bins (in one z bin)
    r0        = dvector(0,N_LENSING-1);
    ds0       = dvector(0,N_LENSING-1);
    err_ds0   = dvector(0,N_LENSING-1);
    res_ds0   = dvector(0,N_LENSING-1);  /* ds pointer: dvector is double, vector is float */
    res_ps0   = dvector(0,N_LENSING-1);
    res_2h0   = dvector(0,N_LENSING-1);
    res_1h0   = dvector(0,N_LENSING-1);
    res_full0 = dvector(0,N_LENSING-1);
    covar0    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r1        = dvector(0,N_LENSING-1);  // data
    ds1       = dvector(0,N_LENSING-1);  // data
    err_ds1   = dvector(0,N_LENSING-1);  // data
    res_ds1   = dvector(0,N_LENSING-1);  // Fit to data, full Delta Sigma
    res_ps1   = dvector(0,N_LENSING-1);  // Fit to data, point source
    res_2h1   = dvector(0,N_LENSING-1);  // Fit to data, two halo
    res_1h1   = dvector(0,N_LENSING-1);  // Fit to data, one halo
    res_full1 = dvector(0,N_LENSING-1);  // Fit to data, full Delta Sigma, using different integration
    covar1    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r2        = dvector(0,N_LENSING-1);
    ds2       = dvector(0,N_LENSING-1);
    err_ds2   = dvector(0,N_LENSING-1);
    res_ds2   = dvector(0,N_LENSING-1);
    res_ps2   = dvector(0,N_LENSING-1);
    res_2h2   = dvector(0,N_LENSING-1);
    res_1h2   = dvector(0,N_LENSING-1);
    res_full2 = dvector(0,N_LENSING-1);
    covar2    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r3        = dvector(0,N_LENSING-1);
    ds3       = dvector(0,N_LENSING-1);
    err_ds3   = dvector(0,N_LENSING-1);
    res_ds3   = dvector(0,N_LENSING-1);
    res_ps3   = dvector(0,N_LENSING-1);
    res_2h3   = dvector(0,N_LENSING-1);
    res_1h3   = dvector(0,N_LENSING-1);
    res_full3 = dvector(0,N_LENSING-1);
    covar3    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r4        = dvector(0,N_LENSING-1);
    ds4       = dvector(0,N_LENSING-1);
    err_ds4   = dvector(0,N_LENSING-1);
    res_ds4   = dvector(0,N_LENSING-1);
    res_ps4   = dvector(0,N_LENSING-1);
    res_2h4   = dvector(0,N_LENSING-1);
    res_1h4   = dvector(0,N_LENSING-1);
    res_full4 = dvector(0,N_LENSING-1);
    covar4    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r5        = dvector(0,N_LENSING-1);
    ds5       = dvector(0,N_LENSING-1);
    err_ds5   = dvector(0,N_LENSING-1);
    res_ds5   = dvector(0,N_LENSING-1);
    res_ps5   = dvector(0,N_LENSING-1);
    res_2h5   = dvector(0,N_LENSING-1);
    res_1h5   = dvector(0,N_LENSING-1);
    res_full5 = dvector(0,N_LENSING-1);
    covar5    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r6        = dvector(0,N_LENSING-1);
    ds6       = dvector(0,N_LENSING-1);
    err_ds6   = dvector(0,N_LENSING-1);
    res_ds6   = dvector(0,N_LENSING-1);
    res_ps6   = dvector(0,N_LENSING-1);
    res_2h6   = dvector(0,N_LENSING-1);
    res_1h6   = dvector(0,N_LENSING-1);
    res_full6 = dvector(0,N_LENSING-1);
    covar6    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    // READ IN THE DATA

    // ---------------------- SM0 ------------------------------------------------
    // Using the /Results_gg/SM_PAPER results here.

    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm0.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm0.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm0.fits.txt");
      }
    }
  
    // 1)mean msun_lens  2)median msun_lens 3)mean z_lens  4)median z_lens  5)plot_radius_kpc[i+1], 
    // 6)radius_kpc[i+1]  7)we1_mean[i+1]  8)e1_num[i]  9)we1_error[i]  10)j_knife_err[i+1]

    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm0_log, &dump, &dump, &dump, &r0[i], &dump, &ds0[i], &dump, &err_ds0[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm0 = pow(10.0,sm0_log);

    // SM0 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm0.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm0.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm0.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar0[i][j]);
      fclose(fp);
      
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar0[i][i] = covar0[i][i]+pow(err_ds0[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar0[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar0[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
    
    // ---------------------- SM1 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm1.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm1.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm1.fits.txt");
      }
    }

    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm1_log, &dump, &dump, &dump, &r1[i], &dump, &ds1[i], &dump, &err_ds1[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm1 = pow(10.0,sm1_log);

    // SM1 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm1.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm1.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm1.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar1[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar1[i][i] = covar1[i][i]+pow(err_ds1[i],2);

       printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar1[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar1[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }

    // ---------------------- SM2 ------------------------------------------------

    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm2.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm2.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm2.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm2_log, &dump, &dump, &dump, &r2[i], &dump, &ds2[i], &dump, &err_ds2[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm2 = pow(10.0,sm2_log);

    // SM2 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm2.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm2.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm2.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar2[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar2[i][i] = covar2[i][i]+pow(err_ds2[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar2[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar2[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
  
    // ---------------------- SM3 ------------------------------------------------
  
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm3.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm3.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm3.fits.txt");
      }
    }
    
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm3_log, &dump, &dump, &dump, &r3[i], &dump, &ds3[i], &dump, &err_ds3[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm3 = pow(10.0,sm3_log);

    // SM3 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm3.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm3.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm3.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar3[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar3[i][i] = covar3[i][i]+pow(err_ds3[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar3[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar3[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
 
    // ---------------------- SM4 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm4.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm4.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm4.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm4_log, &dump, &dump, &dump, &r4[i], &dump, &ds4[i], &dump, &err_ds4[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm4 = pow(10.0,sm4_log);

    // SM4 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm4.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm4.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm4.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar4[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar4[i][i] = covar4[i][i]+pow(err_ds4[i],2);

       printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar4[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar4[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }

    // ---------------------- SM5 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm5.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm5.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm5.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm5_log, &dump, &dump, &dump, &r5[i], &dump, &ds5[i], &dump, &err_ds5[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm5 = pow(10.0,sm5_log);

    // SM5 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm5.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm5.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_blue_sm5.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar5[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar5[i][i] = covar5[i][i]+pow(err_ds5[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar5[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar5[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
  
    // ---------------------- SM6 ------------------------------------------------
  
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm6.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm6.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm6.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm6.fits.txt");
      }
    }

    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      fp = openfile(fname);
      n = filesize(fp);
      for(i=0;i<=N_LENSING-1;i++)
	{
	  //                                              1         2      3      4      5       6      7        8      9            10
	  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm6_log, &dump, &dump, &dump, &r6[i], &dump, &ds6[i], &dump, &err_ds6[i], &dump);  
	  fgets(aa,1000,fp);
	}
      fclose(fp);
      sm6 = pow(10.0,sm6_log);
    }

    // SM6 Lensing covar :
    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

	if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	  sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_blue_sm6.covar");
	if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	  sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_blue_sm6.covar");

	fp = openfile(fname);
	if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	  {
	    fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	    exit(0);
	  }

	for(i=0;i<=N_LENSING-1;++i)
	  for(j=0;j<=N_LENSING-1;++j)
	    fscanf(fp,"%d %d %lf",&i1,&i1,&covar6[i][j]);
	fclose(fp);
    
	// Now add the statistical diagonals before inverting the matrix
	// The covar matrix is the covariance so add err^2
	for(i=0;i<=N_LENSING-1;++i)
	  covar6[i][i] = covar6[i][i]+pow(err_ds0[i],2);

	printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
	tmp=dmatrix(1,N_LENSING,1,1);
	tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
	for(i=1;i<=N_LENSING;++i)
	  for(j=1;j<=N_LENSING;++j)
	    tmp2[i][j]=covar6[i-1][j-1];
	gaussj(tmp2,N_LENSING,tmp,1);
	for(i=1;i<=N_LENSING;++i)
	  for(j=1;j<=N_LENSING;++j)
	    covar6[i-1][j-1]=tmp2[i][j];
	free_dmatrix(tmp,1,N_LENSING,1,1);
	free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
      }
    }

    // ------------------------

    // Data has been read in, set this to 0
    FIRST_LENSING_RUN = 0;
    
    // Test statement
    if(test){
      printf("\n  %f\n",REDSHIFT);   
      for(i=0;i<=N_LENSING-1;i++){
	printf("r, ds, err: %f  %f  %f\n", r1[i],ds1[i],err_ds1[i]);
      }
    } // end of test

    // Test Delta Sigma here (use alog with test_chi2_lensing.pro)
    //test_delta_sigma();

    // TABULATE 2H TO MAKE CODE FASTER
    // !!!!!!!!!!******* NEED TO CHECK THE EFFECTS OF THIS CALCULATION ***** 
    // TABULATE ALL THE BINS HERE (ONLY AT THE BEGINNIN)
    ALEX_2H_BIN_FLAG=0;             // Tabulate the 2h this SM bin
    ALEX_FLAG_2H_TABULATE=1;        // Tabulate the 2h
    RESET_FLAG_2H=1;                // Jeremy flag to reset the 2h
    ALEX_RESET_FLAG_XI=1;           // Alexie flag to retabulate xi
    alex_mmin = sm_lim_1;           // SM bins
    alex_mmax = sm_lim_0;           // SM bins            
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly  
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=1;             // increment here
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1;                         
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=2;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_3;           
    alex_mmax = sm_lim_2;                     
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=3;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_4;           
    alex_mmax = sm_lim_3;               
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);       
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_5;          
    alex_mmax = sm_lim_4;           
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=5;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_6;           
    alex_mmax = sm_lim_5;          
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);       
    dump = tabulate_xi_2h(1.0);

    if(REDSHIFT<Z_CUT3){  // Not for last z bin
      ALEX_2H_BIN_FLAG=6;
      ALEX_FLAG_2H_TABULATE=1;
      RESET_FLAG_2H=1;    
      ALEX_RESET_FLAG_XI=1; 
      alex_mmin = sm_lim_7;          
      alex_mmax = sm_lim_6;                      
      set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      
      dump = tabulate_xi_2h(1.0);
    }
  }// This is the end of : if(FIRST_LENSING_RUN == 1)

  //------------------ END OF READ DATA ------------------------------------------------------

  // Make a pretty plot
  if(LENSING_OUTPUT_FLAG == 3){

    // 2h for EDO
    /*printf("HERE ******* \n");
    ALEX_2H_BIN_FLAG=1;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;  
    alex_mmin = pow(10,  11.4);           
    alex_mmax = pow(10,  11.5); 
    printf("HERE 1 ******* \n");
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);
    dump = tabulate_xi_2h(1.0);
    printf("HERE 2 ******* \n");
    printf("\n >>>>>   >>> REDSHIFT %f \n",REDSHIFT);
    printf("\n >>>>>   >>>mhalo %f \n",ms_to_mhalo(sm_lim_1,wpl.a));
    make_pretty_plot_foredo(sm1);
    exit(0);*/

    // This is to make a pretty plot for paper
    // HIGH MASS (10^11)
    /*ALEX_2H_BIN_FLAG=1;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;          
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    printf("\n pretty_plot_1 : %f \n",sm1_log);
    make_pretty_plot(sm1,1);*/

    // LOW MASS (10^10)
    /*ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;
    alex_mmin = sm_lim_5;           
    alex_mmax = sm_lim_4; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    printf("\n pretty_plot_2 : %f \n",sm4_log);
    make_pretty_plot(sm4,2);*/

    // Sigma for Fabian
    ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;
    //alex_mmin = pow(10, 10.8);  //sm limits here for Fabian           
    //alex_mmax = pow(10, 11.5);
    alex_mmin = pow(10, 11.0);  //sm limits here for ERIC           
    alex_mmax = pow(10, 12.0);
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    sigma_fabian(sm4);

    printf("\n --------- EXIT AFTER MAKE PRETTY PLOT---------------\n");
    exit(0);
  }

  //----------------- Calcualte Delta Sigma --------------------
  // This step takes 0.6s. 0.6*7 = 4.2 (redo timing here ....!)

  // Note: now tabulating 2h
  // Before was doing :
  //RESET_FLAG_2H=1; 
  //dump  = xi_2h_gm(1, alex_mmin, alex_mmax, wpl.a);       // Call this once just to set things up

  // ### CHECK TABULTE 2H FLAGS HERE
  ALEX_2H_BIN_FLAG=0;             // This is for the tabulated 2h
  ALEX_RESET_FLAG_XI=1;
  alex_mmin = sm_lim_1;           // SM bins
  alex_mmax = sm_lim_0;           // SM bins
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
  ds_mod(r0,sm0,res_ds0, res_ps0, res_2h0, res_1h0, res_full0, a);  //#### wpl.a ???? CAN PROBABLY REMOVE THIS a

  ALEX_2H_BIN_FLAG=1;
  ALEX_RESET_FLAG_XI=1;          
  alex_mmin = sm_lim_2;           
  alex_mmax = sm_lim_1;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r1,sm1,res_ds1, res_ps1, res_2h1, res_1h1, res_full1, a);

  ALEX_2H_BIN_FLAG=2;
  ALEX_RESET_FLAG_XI=1;          
  alex_mmin = sm_lim_3;           
  alex_mmax = sm_lim_2;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);  
  ds_mod(r2,sm2,res_ds2, res_ps2, res_2h2, res_1h2, res_full2, a);

  ALEX_2H_BIN_FLAG=3;
  ALEX_RESET_FLAG_XI=1;            
  alex_mmin = sm_lim_4;           
  alex_mmax = sm_lim_3;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r3,sm3,res_ds3, res_ps3, res_2h3, res_1h3, res_full3, a);

  ALEX_2H_BIN_FLAG=4;
  ALEX_RESET_FLAG_XI=1;               
  alex_mmin = sm_lim_5;           
  alex_mmax = sm_lim_4;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);
  ds_mod(r4,sm4,res_ds4, res_ps4, res_2h4, res_1h4, res_full4, a);

  ALEX_2H_BIN_FLAG=5;
  ALEX_RESET_FLAG_XI=1;             
  alex_mmin = sm_lim_6;           
  alex_mmax = sm_lim_5;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r5,sm5,res_ds5, res_ps5, res_2h5, res_1h5, res_full5, a);

  // Not for the last z bin
  if(REDSHIFT<Z_CUT3){
    ALEX_2H_BIN_FLAG=6;  
    ALEX_RESET_FLAG_XI=1;             
    alex_mmin = sm_lim_7;           
    alex_mmax = sm_lim_6;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);   
    ds_mod(r6,sm6,res_ds6, res_ps6, res_2h6, res_1h6, res_full6, a);
  }

  // Write the results to text files to make IDL plots
  if(LENSING_OUTPUT_FLAG){

    // -------- WRITE SM0
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      printf("\n Now writing the lensing results for low z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm0.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      printf("\n Now writing the lensing results for middle z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm0.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      printf("\n Now writing the lensing results for high z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm0.out.txt","w");
    } 

    // Set this to calculate the satellite term
    ALEX_2H_BIN_FLAG=0; 
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_1;           
    alex_mmax = sm_lim_0; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r0[i]);  // Also calculate the 1h satellite term 
      sig_sat_temp = sigma(r0[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r0[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r0[i],res_ds0[i],res_ps0[i],res_2h0[i],res_1h0[i],res_full0[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM1
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm1.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm1.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm1.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=1;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;      
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

   for(i=0;i<=N_LENSING-1;i++){
     ALEX_RESET_FLAG_XI=1; 
     ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
     ds_sat_temp = delta_sigma(r1[i]);  // Also calculate the 1h satellite term
     sig_sat_temp = sigma(r1[i]);
     ALEX_RESET_FLAG_XI=1;  
     ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
     sig_all_temp = sigma(r1[i]);
     fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r1[i],res_ds1[i],res_ps1[i],res_2h1[i],res_1h1[i],res_full1[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM2
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm2.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm2.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm2.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=2;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;      
    alex_mmin = sm_lim_3;           
    alex_mmax = sm_lim_2;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r2[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r2[i]);
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r2[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r2[i],res_ds2[i],res_ps2[i],res_2h2[i],res_1h2[i],res_full2[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM3
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm3.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm3.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm3.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=3;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_4;           
    alex_mmax = sm_lim_3;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r3[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r3[i]);
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r3[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r3[i],res_ds3[i],res_ps3[i],res_2h3[i],res_1h3[i],res_full3[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM4
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm4.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm4.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm4.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=4;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_5;           
    alex_mmax = sm_lim_4; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r4[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r4[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r4[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r4[i],res_ds4[i],res_ps4[i],res_2h4[i],res_1h4[i],res_full4[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM5
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm5.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm5.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_blue_sm5.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=5;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_6;           
    alex_mmax = sm_lim_5;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r5[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r5[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r5[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r5[i],res_ds5[i],res_ps5[i],res_2h5[i],res_1h5[i],res_full5[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }    

    // -------- WRITE SM6
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_blue_sm6.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_blue_sm6.out.txt","w");
    }

    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      ALEX_2H_BIN_FLAG=6;
      ALEX_RESET_FLAG_XI=1;  
      ALEX_XI_OPTION=1; 
      ALEX_PRETTY_PLOT_FLAG=2;     
      alex_mmin = sm_lim_7;           
      alex_mmax = sm_lim_6; 
      set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

      for(i=0;i<=N_LENSING-1;i++){
	ALEX_RESET_FLAG_XI=1; 
	ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
	ds_sat_temp = delta_sigma(r6[i]);  // Also calculate the 1h satellite term
	sig_sat_temp = sigma(r6[i]);
	ALEX_RESET_FLAG_XI=1; 
	ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
	sig_all_temp = sigma(r6[i]);
	fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r6[i],res_ds6[i],res_ps6[i],res_2h6[i],res_1h6[i],res_full6[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
      }
    }
    // printf("\n EXITING IN CHI2_LENSING, 4");
    //exit(0);
  }

  // This computes the chi square
  // Note: data and model are in h72 physical  
  if (LENSING_COVAR){
    // DOUBLE SUM HERE :
    for(i=0;i<=N_LENSING-1;i++){
      for(j=0;j<=N_LENSING-1;j++){
    	chi0=chi0+((res_ds0[i]-ds0[i])*(res_ds0[j]-ds0[j])*covar0[i][j]);
    	chi1=chi1+((res_ds1[i]-ds1[i])*(res_ds1[j]-ds1[j])*covar1[i][j]);
    	chi2=chi2+((res_ds2[i]-ds2[i])*(res_ds2[j]-ds2[j])*covar2[i][j]);
    	chi3=chi3+((res_ds3[i]-ds3[i])*(res_ds3[j]-ds3[j])*covar3[i][j]);
    	chi4=chi4+((res_ds4[i]-ds4[i])*(res_ds4[j]-ds4[j])*covar4[i][j]);
    	chi5=chi5+((res_ds5[i]-ds5[i])*(res_ds5[j]-ds5[j])*covar5[i][j]);
    	if(REDSHIFT<Z_CUT3){
    	  chi6=chi6+((res_ds6[i]-ds6[i])*(res_ds6[j]-ds6[j])*covar6[i][j]);
    	}
      }
    }
  }else{
    // Single sum here :
    for(i=0;i<=N_LENSING-1;i++){
      chi0=chi0+((pow(res_ds0[i]-ds0[i],2))/pow(err_ds0[i],2));
      chi1=chi1+((pow(res_ds1[i]-ds1[i],2))/pow(err_ds1[i],2));
      chi2=chi2+((pow(res_ds2[i]-ds2[i],2))/pow(err_ds2[i],2));
      chi3=chi3+((pow(res_ds3[i]-ds3[i],2))/pow(err_ds3[i],2));
      chi4=chi4+((pow(res_ds4[i]-ds4[i],2))/pow(err_ds4[i],2));
      chi5=chi5+((pow(res_ds5[i]-ds5[i],2))/pow(err_ds5[i],2));
      if(REDSHIFT<Z_CUT3){
	chi6=chi6+((pow(res_ds6[i]-ds6[i],2))/pow(err_ds6[i],2));
      }
    }
  }

  // TINKER wants to take a look
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 0 %e %e\n",r0[i],res_ds0[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 1 %e %e\n",r1[i],res_ds1[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 2 %e %e\n",r2[i],res_ds2[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 3 %e %e\n",r3[i],res_ds3[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 4 %e %e\n",r4[i],res_ds4[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 5 %e %e\n",r5[i],res_ds5[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XX 6 %e %e\n",r6[i],res_ds6[i]);

  // Not for the last z bin
  if(REDSHIFT<Z_CUT3){
    chi_tot = chi0+chi1+chi2+chi3+chi4+chi5+chi6;
  }
  else{
    chi_tot = chi0+chi1+chi2+chi3+chi4+chi5;
  }

  if(test){
    printf("\n");
    printf("chi0 : %f\n", chi0);
    printf("chi1 : %f\n", chi1);
    printf("chi2 : %f\n", chi2);
    printf("chi3 : %f\n", chi3);
    printf("chi4 : %f\n", chi4);
    printf("chi5 : %f\n", chi5);
    printf("chi6 : %f\n", chi6);
    printf("total chi square : %f\n", chi_tot);

    //for(i=0;i<=N_LENSING-1;i++){
    //  printf("\n");
    //  printf("r: %f\n", r4[i]);
    //  printf("ds mod: %f\n", res_ds4[i]);
    //  printf("ds err: %f\n", err_ds4[i]);
    //  printf("ds : %f\n", ds4[i]);
    //}
    exit(0);
  }

  EXCLUSION=4; // Halo exclusion : puts this to galaxy-galaxy cross-correlation for JT

  /*t2=second();
  t1=timediff(t1,t2);
  printf("chi2_lensing : Timing sec: %f\n", t1);*/

  if(ERROR_FLAG == 1){
    printf("\n STOPPED at the bottom of chi2_lensing. There is a previously undetected ERROR_FLAG\n");
    exit(0);
  }

  // printf(" -- NOW EXITING CHI2_LENSING\n");
  return(chi_tot);
}

//-----------------------------------------------------------
//  ***********    READ DATA, CALCULATE CHI2   **************
//  *********************************************************
//-----------------------------------------------------------

// JPL: replace /Users/alexie/Work/
// /home/asleauth/

double chi2_lensing_red(double *a)
{
  int i,j,test,n,i1;
  char fname[1000], aa[1000];                                 // To read in data
  FILE *fp;                                                   // To read in data
  double z,dump;                                              // Redshift and dump vector
  double **tmp,**tmp2;                                        // Dump matrix for inverting covar matrix
  double t1,t2,t3,t4;                                         // Time code
  double ds_sat_temp,sig_sat_temp,sig_cen_temp,sig_all_temp;  // Calcualte 1h_sat for plotting purposes
  double chi0,chi1,chi2,chi3,chi4,chi5,chi6,chi_tot;          // Chi2 for various stellar mass bins
  static int FIRST_LENSING_RUN=1;                             // Set to 1 and put to 0 after data has been read in (!! set to one here to only read in once)
  static float Z_CUT1,Z_CUT2,Z_CUT3,Z_CUT4;                   // Redshift bins
  static double *xx;                                          // dump vector
  static double sm_lim_0,sm_lim_1,sm_lim_2,sm_lim_3,sm_lim_4,sm_lim_5,sm_lim_6,sm_lim_7;  // Stellar mass bins
  static double sm0,sm1,sm2,sm3,sm4,sm5,sm6;                                              // Mean stellar mass in each bin
  static double sm0_log,sm1_log,sm2_log,sm3_log,sm4_log,sm5_log,sm6_log;                  // Log stellar mass
  static double *r0,*ds0,*err_ds0, *res_ds0, *res_ps0, *res_2h0, *res_1h0, *res_full0;    // Data and fits to data
  static double *r1,*ds1,*err_ds1, *res_ds1, *res_ps1, *res_2h1, *res_1h1, *res_full1;
  static double *r2,*ds2,*err_ds2, *res_ds2, *res_ps2, *res_2h2, *res_1h2, *res_full2;
  static double *r3,*ds3,*err_ds3, *res_ds3, *res_ps3, *res_2h3, *res_1h3, *res_full3;
  static double *r4,*ds4,*err_ds4, *res_ds4, *res_ps4, *res_2h4, *res_1h4, *res_full4;
  static double *r5,*ds5,*err_ds5, *res_ds5, *res_ps5, *res_2h5, *res_1h5, *res_full5;
  static double *r6,*ds6,*err_ds6, *res_ds6, *res_ps6, *res_2h6, *res_1h6, *res_full6;
  static double **covar0,**covar1,**covar2,**covar3,**covar4,**covar5,**covar6;           // Covar matrix

  if(DONT_FIT_LENSING)
    return(0.0);

  wpl.reset_inversion = 1;
  BLUE_FLAG = 0;
  EXCLUSION=5;                 // HOD exclusion: this calculates 2h for galaxy-matter cross-correlation 
  test=0;                      // Print some test statements
  N_LENSING =8;                // Number of lensing data points
  ALEX_USE_2H_TABULATE=1;      // Implement the tabulated version of 2h
  ALEX_XI_TEST_OPTION=0;       // Use to test delta_sigma calculation
  
  ALEX_PRETTY_PLOT_FLAG =0;    // (*Don't set here*) Set this later to call 1h_cent and 1h_sat

  //printf(" -- NOW ENTERING CHI2_LENSING\n");
  //t1=second();

  // Use this to make a M-c relation for IDL
  //write_text_halo_concentration();

  // Test Halo Mass Function
  /*  fp=fopen("/Users/alexie/Work/HOD/Pro/dndm_code_0.88.txt","w");
  for(i=0;i<=100;i++){
    fprintf(fp,"%lf %lf \n",10+(5.*i/100.0),log10(dndM_interp(pow(10.0,10+(5.*i/100.0)))));
  }
  fclose(fp);
  exit(0);*/

  if(ERROR_FLAG == 1){
    printf("\n STOPPED at the top of chi2_lensing. There is a previously undetected ERROR_FLAG\n");
    exit(0);
  }

  Z_CUT1 = 0.22;               // Redshift bins
  Z_CUT2 = 0.48; 
  Z_CUT3 = 0.74;
  Z_CUT4 = 1.0;
  
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 7;
  }
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 7;
  }
  if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
    N_LENSING_DATA_BINS = 6;
  }

  // These are the parameters for the final plot
  // This is set in MCMC LENSING
  if(LENSING_OUTPUT_FLAG){
    // Change N_lensing for finner binning (see the files: finebin.txt)
    N_LENSING = 61; 
  }

  chi_tot = 0;
  chi0 = 0;
  chi1 = 0;
  chi2 = 0;
  chi3 = 0;
  chi4 = 0;
  chi5 = 0;
  chi6 = 0;
   
  // READ IN DATA
  // /home/asleauth/Weak_lensing/Results_gg/SM_PAPER
  // NOTE: these distances are: h=0.7, kpc, physical

  if(FIRST_LENSING_RUN == 1){

    // Stellar mass bins
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      //cut1 = [8.70, 9.20, 9.82, 10.30, 10.64, 10.89, 11.12, 12.0]    
      sm_lim_0 = pow(10,  12.0); // hi
      sm_lim_1 = pow(10,  11.12);
      sm_lim_2 = pow(10,  10.89);
      sm_lim_3 = pow(10,  10.64);
      sm_lim_4 = pow(10,  10.30);
      sm_lim_5 = pow(10,  9.82);
      sm_lim_6 = pow(10,  9.20);
      sm_lim_7 = pow(10,  8.70); //low
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      //cut2 = [9.30, 9.80, 10.30, 10.65, 10.88, 11.05, 11.29, 12.0]
      sm_lim_0 = pow(10,  12.0);
      sm_lim_1 = pow(10,  11.29);
      sm_lim_2 = pow(10,  11.05);
      sm_lim_3 = pow(10,  10.88);
      sm_lim_4 = pow(10,  10.65);
      sm_lim_5 = pow(10,  10.30);
      sm_lim_6 = pow(10,  9.80);
      sm_lim_7 = pow(10,  9.30);
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      // NOTE: only  6 bins here
      //cut3 = [9.80, 10.39, 10.74,10.97,11.16,11.35,12.0]
      sm_lim_0 = pow(10,  12.0);
      sm_lim_1 = pow(10,  11.35);
      sm_lim_2 = pow(10,  11.16);
      sm_lim_3 = pow(10,  10.97);
      sm_lim_4 = pow(10,  10.74);
      sm_lim_5 = pow(10,  10.39);
      sm_lim_6 = pow(10,  9.80);
    }

    xx = dvector(0,N_LENSING-1); //dump vector

    // The various lensing stellar mass bins (in one z bin)
    r0        = dvector(0,N_LENSING-1);
    ds0       = dvector(0,N_LENSING-1);
    err_ds0   = dvector(0,N_LENSING-1);
    res_ds0   = dvector(0,N_LENSING-1);  /* ds pointer: dvector is double, vector is float */
    res_ps0   = dvector(0,N_LENSING-1);
    res_2h0   = dvector(0,N_LENSING-1);
    res_1h0   = dvector(0,N_LENSING-1);
    res_full0 = dvector(0,N_LENSING-1);
    covar0    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r1        = dvector(0,N_LENSING-1);  // data
    ds1       = dvector(0,N_LENSING-1);  // data
    err_ds1   = dvector(0,N_LENSING-1);  // data
    res_ds1   = dvector(0,N_LENSING-1);  // Fit to data, full Delta Sigma
    res_ps1   = dvector(0,N_LENSING-1);  // Fit to data, point source
    res_2h1   = dvector(0,N_LENSING-1);  // Fit to data, two halo
    res_1h1   = dvector(0,N_LENSING-1);  // Fit to data, one halo
    res_full1 = dvector(0,N_LENSING-1);  // Fit to data, full Delta Sigma, using different integration
    covar1    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r2        = dvector(0,N_LENSING-1);
    ds2       = dvector(0,N_LENSING-1);
    err_ds2   = dvector(0,N_LENSING-1);
    res_ds2   = dvector(0,N_LENSING-1);
    res_ps2   = dvector(0,N_LENSING-1);
    res_2h2   = dvector(0,N_LENSING-1);
    res_1h2   = dvector(0,N_LENSING-1);
    res_full2 = dvector(0,N_LENSING-1);
    covar2    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r3        = dvector(0,N_LENSING-1);
    ds3       = dvector(0,N_LENSING-1);
    err_ds3   = dvector(0,N_LENSING-1);
    res_ds3   = dvector(0,N_LENSING-1);
    res_ps3   = dvector(0,N_LENSING-1);
    res_2h3   = dvector(0,N_LENSING-1);
    res_1h3   = dvector(0,N_LENSING-1);
    res_full3 = dvector(0,N_LENSING-1);
    covar3    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r4        = dvector(0,N_LENSING-1);
    ds4       = dvector(0,N_LENSING-1);
    err_ds4   = dvector(0,N_LENSING-1);
    res_ds4   = dvector(0,N_LENSING-1);
    res_ps4   = dvector(0,N_LENSING-1);
    res_2h4   = dvector(0,N_LENSING-1);
    res_1h4   = dvector(0,N_LENSING-1);
    res_full4 = dvector(0,N_LENSING-1);
    covar4    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r5        = dvector(0,N_LENSING-1);
    ds5       = dvector(0,N_LENSING-1);
    err_ds5   = dvector(0,N_LENSING-1);
    res_ds5   = dvector(0,N_LENSING-1);
    res_ps5   = dvector(0,N_LENSING-1);
    res_2h5   = dvector(0,N_LENSING-1);
    res_1h5   = dvector(0,N_LENSING-1);
    res_full5 = dvector(0,N_LENSING-1);
    covar5    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    r6        = dvector(0,N_LENSING-1);
    ds6       = dvector(0,N_LENSING-1);
    err_ds6   = dvector(0,N_LENSING-1);
    res_ds6   = dvector(0,N_LENSING-1);
    res_ps6   = dvector(0,N_LENSING-1);
    res_2h6   = dvector(0,N_LENSING-1);
    res_1h6   = dvector(0,N_LENSING-1);
    res_full6 = dvector(0,N_LENSING-1);
    covar6    = dmatrix(0,N_LENSING-1,0,N_LENSING-1);

    // READ IN THE DATA

    // ---------------------- SM0 ------------------------------------------------
    // Using the /Results_gg/SM_PAPER results here.

    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm0.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm0.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm0.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm0.fits.txt");
      }
    }
  
    // 1)mean msun_lens  2)median msun_lens 3)mean z_lens  4)median z_lens  5)plot_radius_kpc[i+1], 
    // 6)radius_kpc[i+1]  7)we1_mean[i+1]  8)e1_num[i]  9)we1_error[i]  10)j_knife_err[i+1]

    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm0_log, &dump, &dump, &dump, &r0[i], &dump, &ds0[i], &dump, &err_ds0[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm0 = pow(10.0,sm0_log);

    // SM0 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm0.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm0.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm0.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar0[i][j]);
      fclose(fp);
      
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar0[i][i] = covar0[i][i]+pow(err_ds0[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar0[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar0[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
    
    // ---------------------- SM1 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm1.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm1.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm1.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm1.fits.txt");
      }
    }

    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm1_log, &dump, &dump, &dump, &r1[i], &dump, &ds1[i], &dump, &err_ds1[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm1 = pow(10.0,sm1_log);

    // SM1 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm1.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm1.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm1.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar1[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar1[i][i] = covar1[i][i]+pow(err_ds1[i],2);

       printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar1[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar1[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }

    // ---------------------- SM2 ------------------------------------------------

    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm2.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm2.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm2.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm2.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm2_log, &dump, &dump, &dump, &r2[i], &dump, &ds2[i], &dump, &err_ds2[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm2 = pow(10.0,sm2_log);

    // SM2 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm2.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm2.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm2.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar2[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar2[i][i] = covar2[i][i]+pow(err_ds2[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar2[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar2[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
  
    // ---------------------- SM3 ------------------------------------------------
  
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm3.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm3.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm3.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm3.fits.txt");
      }
    }
    
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm3_log, &dump, &dump, &dump, &r3[i], &dump, &ds3[i], &dump, &err_ds3[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm3 = pow(10.0,sm3_log);

    // SM3 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm3.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm3.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm3.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar3[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar3[i][i] = covar3[i][i]+pow(err_ds3[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar3[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar3[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
 
    // ---------------------- SM4 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm4.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm4.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm4.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm4.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm4_log, &dump, &dump, &dump, &r4[i], &dump, &ds4[i], &dump, &err_ds4[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm4 = pow(10.0,sm4_log);

    // SM4 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm4.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm4.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm4.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar4[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar4[i][i] = covar4[i][i]+pow(err_ds4[i],2);

       printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar4[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar4[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }

    // ---------------------- SM5 ------------------------------------------------
    
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm5.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm5.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      if(LENSING_OUTPUT_FLAG){
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm5.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm5.fits.txt");
      }
    }
  
    fp = openfile(fname);
    n = filesize(fp);
    for(i=0;i<=N_LENSING-1;i++)
      {
	//                                              1         2      3      4      5       6      7        8      9            10
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm5_log, &dump, &dump, &dump, &r5[i], &dump, &ds5[i], &dump, &err_ds5[i], &dump);  
	fgets(aa,1000,fp);
      }
    fclose(fp);
    sm5 = pow(10.0,sm5_log);

    // SM5 Lensing covar :
    if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

      if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm5.covar");
      if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm5.covar");
      if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4)
	sprintf(fname,"/Users/alexie/Work/HOD/Covar/z3_red_sm5.covar");

      fp = openfile(fname);
      if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	{
	  fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	  exit(0);
	}

      for(i=0;i<=N_LENSING-1;++i)
	for(j=0;j<=N_LENSING-1;++j)
	  fscanf(fp,"%d %d %lf",&i1,&i1,&covar5[i][j]);
      fclose(fp);
    
      // Now add the statistical diagonals before inverting the matrix
      // The covar matrix is the covariance so add err^2
      for(i=0;i<=N_LENSING-1;++i)
	covar5[i][i] = covar5[i][i]+pow(err_ds5[i],2);

      printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
      tmp=dmatrix(1,N_LENSING,1,1);
      tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  tmp2[i][j]=covar5[i-1][j-1];
      gaussj(tmp2,N_LENSING,tmp,1);
      for(i=1;i<=N_LENSING;++i)
	for(j=1;j<=N_LENSING;++j)
	  covar5[i-1][j-1]=tmp2[i][j];
      free_dmatrix(tmp,1,N_LENSING,1,1);
      free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
    }
  
    // ---------------------- SM6 ------------------------------------------------
  
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      if(LENSING_OUTPUT_FLAG){      
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm6.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm6.fits.txt");
      }
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      if(LENSING_OUTPUT_FLAG){ 
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm6.fits.finebin.txt");
      }
      else{
	sprintf(fname,"/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm6.fits.txt");
      }
    }

    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      fp = openfile(fname);
      n = filesize(fp);
      for(i=0;i<=N_LENSING-1;i++)
	{
	  //                                              1         2      3      4      5       6      7        8      9            10
	  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&sm6_log, &dump, &dump, &dump, &r6[i], &dump, &ds6[i], &dump, &err_ds6[i], &dump);  
	  fgets(aa,1000,fp);
	}
      fclose(fp);
      sm6 = pow(10.0,sm6_log);
    }

    // SM6 Lensing covar :
    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      if(!LENSING_OUTPUT_FLAG && LENSING_COVAR){ 

	if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2)
	  sprintf(fname,"/Users/alexie/Work/HOD/Covar/z1_red_sm6.covar");
	if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3)
	  sprintf(fname,"/Users/alexie/Work/HOD/Covar/z2_red_sm6.covar");

	fp = openfile(fname);
	if(filesize(fp)!=(N_LENSING)*(N_LENSING))
	  {
	    fprintf(stderr,"FILESIZE mismatch: cij and data not match in CHI2 LENSING: %d %d\n",(N_LENSING)*(N_LENSING),filesize(fp));
	    exit(0);
	  }

	for(i=0;i<=N_LENSING-1;++i)
	  for(j=0;j<=N_LENSING-1;++j)
	    fscanf(fp,"%d %d %lf",&i1,&i1,&covar6[i][j]);
	fclose(fp);
    
	// Now add the statistical diagonals before inverting the matrix
	// The covar matrix is the covariance so add err^2
	for(i=0;i<=N_LENSING-1;++i)
	  covar6[i][i] = covar6[i][i]+pow(err_ds0[i],2);

	printf("INVERTING LENSING LENSING_COVARIANCE MATRIX\n");
	tmp=dmatrix(1,N_LENSING,1,1);
	tmp2=dmatrix(1,N_LENSING,1,N_LENSING);
	for(i=1;i<=N_LENSING;++i)
	  for(j=1;j<=N_LENSING;++j)
	    tmp2[i][j]=covar6[i-1][j-1];
	gaussj(tmp2,N_LENSING,tmp,1);
	for(i=1;i<=N_LENSING;++i)
	  for(j=1;j<=N_LENSING;++j)
	    covar6[i-1][j-1]=tmp2[i][j];
	free_dmatrix(tmp,1,N_LENSING,1,1);
	free_dmatrix(tmp2,1,N_LENSING,1,N_LENSING);
      }
    }

    // ------------------------

    // Data has been read in, set this to 0
    FIRST_LENSING_RUN = 0;
    
    // Test statement
    if(test){
      printf("\n  %f\n",REDSHIFT);   
      for(i=0;i<=N_LENSING-1;i++){
	printf("r, ds, err: %f  %f  %f\n", r1[i],ds1[i],err_ds1[i]);
      }
    } // end of test

    // Test Delta Sigma here (use alog with test_chi2_lensing.pro)
    //test_delta_sigma();

    // TABULATE 2H TO MAKE CODE FASTER
    // !!!!!!!!!!******* NEED TO CHECK THE EFFECTS OF THIS CALCULATION ***** 
    // TABULATE ALL THE BINS HERE (ONLY AT THE BEGINNIN)
    ALEX_2H_BIN_FLAG=0;             // Tabulate the 2h this SM bin
    ALEX_FLAG_2H_TABULATE=1;        // Tabulate the 2h
    RESET_FLAG_2H=1;                // Jeremy flag to reset the 2h
    ALEX_RESET_FLAG_XI=1;           // Alexie flag to retabulate xi
    alex_mmin = sm_lim_1;           // SM bins
    alex_mmax = sm_lim_0;           // SM bins            
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly  
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=1;             // increment here
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1;                         
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=2;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_3;           
    alex_mmax = sm_lim_2;                     
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=3;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_4;           
    alex_mmax = sm_lim_3;               
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);       
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_5;          
    alex_mmax = sm_lim_4;           
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     
    dump = tabulate_xi_2h(1.0);

    ALEX_2H_BIN_FLAG=5;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1; 
    alex_mmin = sm_lim_6;           
    alex_mmax = sm_lim_5;          
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);       
    dump = tabulate_xi_2h(1.0);

    if(REDSHIFT<Z_CUT3){  // Not for last z bin
      ALEX_2H_BIN_FLAG=6;
      ALEX_FLAG_2H_TABULATE=1;
      RESET_FLAG_2H=1;    
      ALEX_RESET_FLAG_XI=1; 
      alex_mmin = sm_lim_7;          
      alex_mmax = sm_lim_6;                      
      set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      
      dump = tabulate_xi_2h(1.0);
    }
  }// This is the end of : if(FIRST_LENSING_RUN == 1)

  //------------------ END OF READ DATA ------------------------------------------------------

  // Make a pretty plot
  if(LENSING_OUTPUT_FLAG == 3){

    // 2h for EDO
    /*printf("HERE ******* \n");
    ALEX_2H_BIN_FLAG=1;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;  
    alex_mmin = pow(10,  11.4);           
    alex_mmax = pow(10,  11.5); 
    printf("HERE 1 ******* \n");
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);
    dump = tabulate_xi_2h(1.0);
    printf("HERE 2 ******* \n");
    printf("\n >>>>>   >>> REDSHIFT %f \n",REDSHIFT);
    printf("\n >>>>>   >>>mhalo %f \n",ms_to_mhalo(sm_lim_1,wpl.a));
    make_pretty_plot_foredo(sm1);
    exit(0);*/

    // This is to make a pretty plot for paper
    // HIGH MASS (10^11)
    /*ALEX_2H_BIN_FLAG=1;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;          
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    printf("\n pretty_plot_1 : %f \n",sm1_log);
    make_pretty_plot(sm1,1);*/

    // LOW MASS (10^10)
    /*ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;
    alex_mmin = sm_lim_5;           
    alex_mmax = sm_lim_4; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    printf("\n pretty_plot_2 : %f \n",sm4_log);
    make_pretty_plot(sm4,2);*/

    // Sigma for Fabian
    ALEX_2H_BIN_FLAG=4;
    ALEX_FLAG_2H_TABULATE=1;
    RESET_FLAG_2H=1;    
    ALEX_RESET_FLAG_XI=1;
    //alex_mmin = pow(10, 10.8);  //sm limits here for Fabian           
    //alex_mmax = pow(10, 11.5);
    alex_mmin = pow(10, 11.0);  //sm limits here for ERIC           
    alex_mmax = pow(10, 12.0);
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
    sigma_fabian(sm4);

    printf("\n --------- EXIT AFTER MAKE PRETTY PLOT---------------\n");
    exit(0);
  }

  //----------------- Calcualte Delta Sigma --------------------
  // This step takes 0.6s. 0.6*7 = 4.2 (redo timing here ....!)

  // Note: now tabulating 2h
  // Before was doing :
  //RESET_FLAG_2H=1; 
  //dump  = xi_2h_gm(1, alex_mmin, alex_mmax, wpl.a);       // Call this once just to set things up

  // ### CHECK TABULTE 2H FLAGS HERE
  ALEX_2H_BIN_FLAG=0;             // This is for the tabulated 2h
  ALEX_RESET_FLAG_XI=1;
  alex_mmin = sm_lim_1;           // SM bins
  alex_mmax = sm_lim_0;           // SM bins
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly
  ds_mod(r0,sm0,res_ds0, res_ps0, res_2h0, res_1h0, res_full0, a);  //#### wpl.a ???? CAN PROBABLY REMOVE THIS a

  ALEX_2H_BIN_FLAG=1;
  ALEX_RESET_FLAG_XI=1;          
  alex_mmin = sm_lim_2;           
  alex_mmax = sm_lim_1;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r1,sm1,res_ds1, res_ps1, res_2h1, res_1h1, res_full1, a);

  ALEX_2H_BIN_FLAG=2;
  ALEX_RESET_FLAG_XI=1;          
  alex_mmin = sm_lim_3;           
  alex_mmax = sm_lim_2;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);  
  ds_mod(r2,sm2,res_ds2, res_ps2, res_2h2, res_1h2, res_full2, a);

  ALEX_2H_BIN_FLAG=3;
  ALEX_RESET_FLAG_XI=1;            
  alex_mmin = sm_lim_4;           
  alex_mmax = sm_lim_3;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r3,sm3,res_ds3, res_ps3, res_2h3, res_1h3, res_full3, a);

  ALEX_2H_BIN_FLAG=4;
  ALEX_RESET_FLAG_XI=1;               
  alex_mmin = sm_lim_5;           
  alex_mmax = sm_lim_4;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);
  ds_mod(r4,sm4,res_ds4, res_ps4, res_2h4, res_1h4, res_full4, a);

  ALEX_2H_BIN_FLAG=5;
  ALEX_RESET_FLAG_XI=1;             
  alex_mmin = sm_lim_6;           
  alex_mmax = sm_lim_5;
  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 
  ds_mod(r5,sm5,res_ds5, res_ps5, res_2h5, res_1h5, res_full5, a);

  // Not for the last z bin
  if(REDSHIFT<Z_CUT3){
    ALEX_2H_BIN_FLAG=6;  
    ALEX_RESET_FLAG_XI=1;             
    alex_mmin = sm_lim_7;           
    alex_mmax = sm_lim_6;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);   
    ds_mod(r6,sm6,res_ds6, res_ps6, res_2h6, res_1h6, res_full6, a);
  }

  // Write the results to text files to make IDL plots
  if(LENSING_OUTPUT_FLAG){

    // -------- WRITE SM0
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      printf("\n Now writing the lensing results for low z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm0.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      printf("\n Now writing the lensing results for middle z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm0.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      printf("\n Now writing the lensing results for high z bin ....");
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm0.out.txt","w");
    } 

    // Set this to calculate the satellite term
    ALEX_2H_BIN_FLAG=0; 
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_1;           
    alex_mmax = sm_lim_0; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a); 

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r0[i]);  // Also calculate the 1h satellite term 
      sig_sat_temp = sigma(r0[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r0[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r0[i],res_ds0[i],res_ps0[i],res_2h0[i],res_1h0[i],res_full0[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM1
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm1.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm1.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm1.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=1;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;      
    alex_mmin = sm_lim_2;           
    alex_mmax = sm_lim_1; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

   for(i=0;i<=N_LENSING-1;i++){
     ALEX_RESET_FLAG_XI=1; 
     ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
     ds_sat_temp = delta_sigma(r1[i]);  // Also calculate the 1h satellite term
     sig_sat_temp = sigma(r1[i]);
     ALEX_RESET_FLAG_XI=1;  
     ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
     sig_all_temp = sigma(r1[i]);
     fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r1[i],res_ds1[i],res_ps1[i],res_2h1[i],res_1h1[i],res_full1[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM2
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm2.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm2.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm2.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=2;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;      
    alex_mmin = sm_lim_3;           
    alex_mmax = sm_lim_2;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r2[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r2[i]);
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r2[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r2[i],res_ds2[i],res_ps2[i],res_2h2[i],res_1h2[i],res_full2[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM3
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm3.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm3.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm3.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=3;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_4;           
    alex_mmax = sm_lim_3;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r3[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r3[i]);
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r3[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r3[i],res_ds3[i],res_ps3[i],res_2h3[i],res_1h3[i],res_full3[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM4
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm4.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm4.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm4.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=4;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_5;           
    alex_mmax = sm_lim_4; 
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r4[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r4[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r4[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r4[i],res_ds4[i],res_ps4[i],res_2h4[i],res_1h4[i],res_full4[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }

    // -------- WRITE SM5
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm5.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm5.out.txt","w");
    }
    if(REDSHIFT>Z_CUT3 && REDSHIFT<Z_CUT4){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z3_red_sm5.out.txt","w");
    }
    ALEX_2H_BIN_FLAG=5;
    ALEX_RESET_FLAG_XI=1;  
    ALEX_XI_OPTION=1; 
    ALEX_PRETTY_PLOT_FLAG=2;     
    alex_mmin = sm_lim_6;           
    alex_mmax = sm_lim_5;
    set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      

    for(i=0;i<=N_LENSING-1;i++){
      ALEX_RESET_FLAG_XI=1; 
      ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
      ds_sat_temp = delta_sigma(r5[i]);  // Also calculate the 1h satellite term
      sig_sat_temp = sigma(r5[i]);
      ALEX_RESET_FLAG_XI=1;  
      ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
      sig_all_temp = sigma(r5[i]);
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r5[i],res_ds5[i],res_ps5[i],res_2h5[i],res_1h5[i],res_full5[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
    }    

    // -------- WRITE SM6
    if(REDSHIFT>Z_CUT1 && REDSHIFT<Z_CUT2){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z1_red_sm6.out.txt","w");
    }
    if(REDSHIFT>Z_CUT2 && REDSHIFT<Z_CUT3){
      fp=fopen("/Users/tinker/cosmo/COSMOS/QUENCHED_LENSING/z2_red_sm6.out.txt","w");
    }

    // Not for the last z bin
    if(REDSHIFT<Z_CUT3){
      ALEX_2H_BIN_FLAG=6;
      ALEX_RESET_FLAG_XI=1;  
      ALEX_XI_OPTION=1; 
      ALEX_PRETTY_PLOT_FLAG=2;     
      alex_mmin = sm_lim_7;           
      alex_mmax = sm_lim_6; 
      set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);     

      for(i=0;i<=N_LENSING-1;i++){
	ALEX_RESET_FLAG_XI=1; 
	ALEX_PRETTY_PLOT_FLAG = 2; // Only sat
	ds_sat_temp = delta_sigma(r6[i]);  // Also calculate the 1h satellite term
	sig_sat_temp = sigma(r6[i]);
	ALEX_RESET_FLAG_XI=1; 
	ALEX_PRETTY_PLOT_FLAG = 0; // Cen and Sat together      
	sig_all_temp = sigma(r6[i]);
	fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",r6[i],res_ds6[i],res_ps6[i],res_2h6[i],res_1h6[i],res_full6[i],ds_sat_temp,sig_sat_temp,sig_all_temp);
      }
    }
    //    printf("\n EXITING IN CHI2_LENSING, 4");
    //exit(0);
  }

  // This computes the chi square
  // Note: data and model are in h72 physical  
  if (LENSING_COVAR){
    // DOUBLE SUM HERE :
    for(i=0;i<=N_LENSING-1;i++){
      for(j=0;j<=N_LENSING-1;j++){
    	chi0=chi0+((res_ds0[i]-ds0[i])*(res_ds0[j]-ds0[j])*covar0[i][j]);
    	chi1=chi1+((res_ds1[i]-ds1[i])*(res_ds1[j]-ds1[j])*covar1[i][j]);
    	chi2=chi2+((res_ds2[i]-ds2[i])*(res_ds2[j]-ds2[j])*covar2[i][j]);
    	chi3=chi3+((res_ds3[i]-ds3[i])*(res_ds3[j]-ds3[j])*covar3[i][j]);
    	chi4=chi4+((res_ds4[i]-ds4[i])*(res_ds4[j]-ds4[j])*covar4[i][j]);
    	chi5=chi5+((res_ds5[i]-ds5[i])*(res_ds5[j]-ds5[j])*covar5[i][j]);
    	if(REDSHIFT<Z_CUT3){
    	  chi6=chi6+((res_ds6[i]-ds6[i])*(res_ds6[j]-ds6[j])*covar6[i][j]);
    	}
      }
    }
  }else{
    // Single sum here :
    for(i=0;i<=N_LENSING-1;i++){
      chi0=chi0+((pow(res_ds0[i]-ds0[i],2))/pow(err_ds0[i],2));
      chi1=chi1+((pow(res_ds1[i]-ds1[i],2))/pow(err_ds1[i],2));
      chi2=chi2+((pow(res_ds2[i]-ds2[i],2))/pow(err_ds2[i],2));
      chi3=chi3+((pow(res_ds3[i]-ds3[i],2))/pow(err_ds3[i],2));
      chi4=chi4+((pow(res_ds4[i]-ds4[i],2))/pow(err_ds4[i],2));
      chi5=chi5+((pow(res_ds5[i]-ds5[i],2))/pow(err_ds5[i],2));
      if(REDSHIFT<Z_CUT3){
	chi6=chi6+((pow(res_ds6[i]-ds6[i],2))/pow(err_ds6[i],2));
      }
    }
  }

  // TINKER wants to take a look
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 0 %e %e\n",r0[i],res_ds0[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 1 %e %e\n",r1[i],res_ds1[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 2 %e %e\n",r2[i],res_ds2[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 3 %e %e\n",r3[i],res_ds3[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 4 %e %e\n",r4[i],res_ds4[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 5 %e %e\n",r5[i],res_ds5[i]);
  for(i=1;i<=N_LENSING-1;++i)
    printf("XR 6 %e %e\n",r6[i],res_ds6[i]);


  // TINKER SETS the SM0 RED CHI^2 to ZERO
  //----------------------------------------
  chi0 = 0;

  // Not for the last z bin
  if(REDSHIFT<Z_CUT3){
    chi_tot = chi0+chi1+chi2+chi3+chi4+chi5+chi6;
  }
  else{
    chi_tot = chi0+chi1+chi2+chi3+chi4+chi5;
  }

  if(test){
    printf("\n");
    printf("chi0 : %f\n", chi0);
    printf("chi1 : %f\n", chi1);
    printf("chi2 : %f\n", chi2);
    printf("chi3 : %f\n", chi3);
    printf("chi4 : %f\n", chi4);
    printf("chi5 : %f\n", chi5);
    printf("chi6 : %f\n", chi6);
    printf("total chi square : %f\n", chi_tot);

    //for(i=0;i<=N_LENSING-1;i++){
    //  printf("\n");
    //  printf("r: %f\n", r4[i]);
    //  printf("ds mod: %f\n", res_ds4[i]);
    //  printf("ds err: %f\n", err_ds4[i]);
    //  printf("ds : %f\n", ds4[i]);
    //}
    exit(0);
  }

  EXCLUSION=4; // Halo exclusion : puts this to galaxy-galaxy cross-correlation for JT

  /*t2=second();
  t1=timediff(t1,t2);
  printf("chi2_lensing : Timing sec: %f\n", t1);*/

  if(ERROR_FLAG == 1){
    printf("\n STOPPED at the bottom of chi2_lensing. There is a previously undetected ERROR_FLAG\n");
    exit(0);
  }

  // printf(" -- NOW EXITING CHI2_LENSING\n");
  return(chi_tot);
}


//-----------------------------------------------------------
//  ***********   CALCULATE XI   ****************************
//  *********************************************************
//-----------------------------------------------------------

//-----------------------------------------------------------
// Tabulate xi (for the first part of the integral)
//-----------------------------------------------------------

double tabulate_xi(double r){
  static int flag=0;        // static means that the function remembers the values
  static double *x,*y,*y2;  // r,xi, and derivative
  int i,n=200;              // number of tabulated points (decrease this to speed up code). 500 is too large. 100 looks good.
  double a;

  if(!flag || ALEX_RESET_FLAG_XI)
    {
      if(!flag) // only allocate memory once for all the bins
	{
	  // from 1 to n size vectors
	  // x and y in linear units, x in linear Mpc
	  x  = dvector(1,n);  // log spaced tabulations
	  y  = dvector(1,n);
	  y2 = dvector(1,n);
	}
      flag=1;

      ALEX_RESET_FLAG_XI=0;
      calc_xi(x,y,n);
      // y2 is the derivative
      // can use in either lin or log space
      // this function sets up the derivatives
      spline(x,y,n,2.0E+30,2.0E+30,y2);  // spine in linear space (but bins are logarithmic)
    }

  // Interpolate beyond limits
  // Manually set to the innermost and outmost value
  //if(r>x[n]) return(y[n]);
  //if(r<x[1]) return(y[1]);

  // Interpolate within limits
  splint(x,y,y2,n,r,&a);
  return(a);
}

//-----------------------------------------------------------
// Tabulate xi_2h_gm (to make the code go faster)
// Should do this separately for all SM bins but for now
// Just do for one SM bin and apply to all others
// COSMOS results are not very sensitive to the 2h
//-----------------------------------------------------------

double tabulate_xi_2h(double r) 
{
  static int flag=0;            // static means that the function remembers the values
  static double *x0,*y0,*y2_0;  // r,xi, and derivative
  static double *x1,*y1,*y2_1;
  static double *x2,*y2,*y2_2;
  static double *x3,*y3,*y2_3;
  static double *x4,*y4,*y2_4;
  static double *x5,*y5,*y2_5;
  static double *x6,*y6,*y2_6;

  int i,n=200;              // number of tabulated points
  double a;
  double rhi,rlo,dr;

  a=0;  // Result 

  if(!flag || ALEX_FLAG_2H_TABULATE)
    {
      if(!flag) // only allocate memory once (FOR ALL THE BINS)
	{
	  // from 1 to n size vectors
	  // x and y in linear units, x in linear Mpc
	  x0   = dvector(1,n);  
	  y0   = dvector(1,n);
	  y2_0 = dvector(1,n);

	  x1   = dvector(1,n);  
	  y1   = dvector(1,n);
	  y2_1 = dvector(1,n);

	  x2   = dvector(1,n);  
	  y2   = dvector(1,n);
	  y2_2 = dvector(1,n);

	  x3   = dvector(1,n);  
	  y3   = dvector(1,n);
	  y2_3 = dvector(1,n);

	  x4   = dvector(1,n);  
	  y4   = dvector(1,n);
	  y2_4 = dvector(1,n);

	  x5   = dvector(1,n);  
	  y5   = dvector(1,n);
	  y2_5 = dvector(1,n);

	  x6   = dvector(1,n);  
	  y6   = dvector(1,n);
	  y2_6 = dvector(1,n);
	}
      rlo=log(1e-4);        // Minimum value
      rhi=log(70);          // Maximum value
      dr=(rhi-rlo)/(n-1);   // Bin size
      flag=1;
      ALEX_FLAG_2H_TABULATE=0;

      // Calculate the 2h term for each of the 7 bins

      //SM0
      if(ALEX_2H_BIN_FLAG==0){
	for(i=1;i<=n;++i)
	  {
	    x0[i]=exp((i-1)*dr + rlo); // x=R is in non log values, Mpc
	    y0[i]=xi_2h_gm(x0[i], alex_mmin, alex_mmax, wpl.a); // y=xi Return linear values
	  }
	// y2 is the derivative
	// can use in either lin or log space
	// this function sets up the derivatives
	spline(x0,y0,n,2.0E+30,2.0E+30,y2_0);  // spine in linear space
      }
	
      //SM1
      if(ALEX_2H_BIN_FLAG==1){
	for(i=1;i<=n;++i)
	  {
	    x1[i]=exp((i-1)*dr + rlo); 
	    y1[i]=xi_2h_gm(x1[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x1,y1,n,2.0E+30,2.0E+30,y2_1); 
      }

      //SM2
      if(ALEX_2H_BIN_FLAG==2){
	for(i=1;i<=n;++i)
	  {
	    x2[i]=exp((i-1)*dr + rlo); 
	    y2[i]=xi_2h_gm(x2[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x2,y2,n,2.0E+30,2.0E+30,y2_2); 
      }

      //SM3
      if(ALEX_2H_BIN_FLAG==3){
	for(i=1;i<=n;++i)
	  {
	    x3[i]=exp((i-1)*dr + rlo); 
	    y3[i]=xi_2h_gm(x3[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x3,y3,n,2.0E+30,2.0E+30,y2_3); 
      }

      //SM4
      if(ALEX_2H_BIN_FLAG==4){
	for(i=1;i<=n;++i)
	  {
	    x4[i]=exp((i-1)*dr + rlo); 
	    y4[i]=xi_2h_gm(x4[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x4,y4,n,2.0E+30,2.0E+30,y2_4); 
      }

      //SM5
      if(ALEX_2H_BIN_FLAG==5){
	for(i=1;i<=n;++i)
	  {
	    x5[i]=exp((i-1)*dr + rlo); 
	    y5[i]=xi_2h_gm(x5[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x5,y5,n,2.0E+30,2.0E+30,y2_5); 
      }

      //SM6
      if(ALEX_2H_BIN_FLAG==6){
	for(i=1;i<=n;++i)
	  {
	    x6[i]=exp((i-1)*dr + rlo); 
	    y6[i]=xi_2h_gm(x6[i], alex_mmin, alex_mmax, wpl.a); 
	  }
	spline(x6,y6,n,2.0E+30,2.0E+30,y2_6); 
      }
      
    }

  // Interpolate within limits

  // SM0
  if(ALEX_2H_BIN_FLAG==0)
    splint(x0,y0,y2_0,n,r,&a);
  
  // SM1
  if(ALEX_2H_BIN_FLAG==1)
    splint(x1,y1,y2_1,n,r,&a);

  // SM2
  if(ALEX_2H_BIN_FLAG==2)
    splint(x2,y2,y2_2,n,r,&a);

  // SM3
  if(ALEX_2H_BIN_FLAG==3)
    splint(x3,y3,y2_3,n,r,&a);

  // SM4
  if(ALEX_2H_BIN_FLAG==4)
    splint(x4,y4,y2_4,n,r,&a);

  // SM5
  if(ALEX_2H_BIN_FLAG==5)
    splint(x5,y5,y2_5,n,r,&a);

  // SM6
  if(ALEX_2H_BIN_FLAG==6)
    splint(x6,y6,y2_6,n,r,&a);

  return(a);
}

//-----------------------------------------------------------
// Calculate xi logarithmically spaced in r
// Returns xi in linear units and r in linear Mpc
//-----------------------------------------------------------

// Returns 1+xi: in units of co h^-1 MPC
void calc_xi(double *r, double *xi, int n)
{
  double rhi,rlo,dr;
  int i;
  FILE *fp;

  // rlo=log(1.0e-4);        // Minimum value (co h^-1 MPC) 
  rlo = -9.210305;// TINKER
  //  rhi=log(70.0);          // Maximum value (co h^-1 MPC)
  rhi = 4.248495e+00;
  dr=(rhi-rlo)/(n-1);   // Logarithmic bin size

  /*
  fmuh(log(0.001));
  fmuh(log(70.0));
  fmuh(rlo);
  fmuh(rhi);
  fmuh(rlo);
  fmuh(rlo);
  fmuh(rhi);
  */

  for(i=1;i<=n;++i)
    {
      r[i]=exp((i-1)*dr + rlo); // R is in non log values, Mpc

      // One Halo xi (calculate 1+x1) 
      // This output is already 1+xi
      if(ALEX_XI_OPTION == 1){
	if(ALEX_XI_TEST_OPTION ==1){
	  xi[i]= xi_1h_simple(r[i]);  // This calculates a simple NFW
	} else{
	  xi[i]= xi_1h(r[i]);         // Return linear values. 
	}
      }

      // Two halo xi (calculate 1+x1)
      if(ALEX_XI_OPTION == 2){
	if(ALEX_USE_2H_TABULATE ==1){
	  xi[i]=1+tabulate_xi_2h(r[i]); // Return linear values
	}
	else{
	  xi[i]=1+xi_2h_gm(r[i], alex_mmin, alex_mmax, wpl.a); // Return linear values //#### check this function
	}
      }

      // Full model: integrate both xi together //## +1 0r +2 here ???
      if(ALEX_XI_OPTION == 3){
	if(ALEX_USE_2H_TABULATE ==1){
	  xi[i]= 1 + xi_1h(r[i]) + tabulate_xi_2h(r[i]);
	}
	else{
	  xi[i]= 1 + xi_1h(r[i]) + xi_2h_gm(r[i], alex_mmin, alex_mmax, wpl.a);
	}
      }

      //      fprintf(stdout,"calc_xi> %d %e %e %d %d %e %e\n",i,r[i],xi[i],ALEX_XI_OPTION,n,dr,rlo);
      //fflush(stdout);
    }
  //exit(0);
}

//-----------------------------------------------------------
// Calculate xi for the One Halo (see Eq 6 Yoo et al.)
// this is called by <calc_xi>
//-----------------------------------------------------------

// Input R is in co H=100
// Output xi is in co H=100 
// Output is ALREADY 1+xi
double xi_1h(double r)
{
  int i;
  double res,norm,ng,mlo,mmin;
  double temp, m,r_vir,x,fsat,fcen,conc,n_sat,n_cen,n_halo,rho_mean;
  double testr,delta,test_val;
  double test_int;

  // NB! setting an upper limit on 1 halo integral.
  // Only cent
  if(ALEX_PRETTY_PLOT_FLAG == 1){
    if(r>1)return 0;
  } else{
    if(r>10)return 0; // Used to be 8 !!!!  ****
  }

  alex_sat_r=r;                          // Calculate at this radius, global variable (this is a 3d dis), dont confuse with alex_trans_R
  ng = GALAXY_DENSITY;                   // ng, galaxy density
  norm=1.0/(4*PI*pow(alex_sat_r,2)*ng);  // normalization in front of equation, Eq 6, 1/(4*pi*r^2*ng)

  // Only integrate from m low
  // Because F(x) goes to zero for x>1
  mlo = HOD.M_low;

  // Make a test integral run
  // To avoid zero integrals
  // Qromo has trouble integrating when the result is zero
  test_int = 0.0;
  delta = (log(HOD.M_max)-log(mlo))/49.0;
  if(isnan(delta))
    {
      delta = (log(HOD.M_max)-log(mlo))/49.0;
      testr = (log(HOD.M_max)-log(mlo))/49.0;
      //printf("isnan> %e %e\n",delta,testr);
      if(!isnan(testr))
	delta = testr;
    }

  for(i=0;i<=49;i++){
    testr  = ((i*delta)+log(mlo));               // Test radius in log units
    //printf("testr> %e %e %e %e %e\n",delta,mlo,testr,HOD.M_max,(log(HOD.M_max)-log(mlo))/49.0);
    test_val = funct_1h(testr);
    test_int = test_int + fabs(test_val);
    // printf("\n  >>> %f %f %f :",testr,test_val,test_int);
  }
  res = 0.0;
  
  if((mlo<HOD.M_max) && (test_int > 1e-8)){ 
    //printf("xi_1h> %e %e\n",log(mlo),log(HOD.M_max));
    //fflush(stdout);
    //log(mlo); log(HOD.M_max);
    //printf("");
    //res = norm*qtrap(funct_1h,log(mlo),log(HOD.M_max),1.0E-4);
    res = norm*qromo(funct_1h,log(mlo),log(HOD.M_max),alex1midpnt);  // integrate in ln space
  }

  //  printf("xi_1h> %e %e %e %e\n",r,mlo,HOD.M_max,res);

  /*  if(mlo>=HOD.M_max){
    res = 0;  
    }*/

  if(ERROR_FLAG == 1){
    printf("\n STOPPED IN xi_1h \n",res);
    printf("\n TEST_INT %f \n",test_int);
    printf(" --------- \n"); 
    printf(" mlo in log10 %f \n",log(mlo)/log(10));  
    printf(" M_max in log10 %f \n",log(HOD.M_max)/log(10));   
    printf(" res %f \n",res);
    printf("norm %f\n",norm);
    printf("alex_sat_r %f\n",alex_sat_r);
    printf("r %f\n",r);
    printf("ng %f\n",ng);
    printf("funct_1h evaluated at mlo %e \n", funct_1h(log(mlo)));
    printf("funct_1h evaluated at mmax %e \n", funct_1h(log(HOD.M_max)));

    printf("\n a[1] %f\n",wpl.a[1]);
    printf("a[2] %f\n",wpl.a[2]);
    printf("a[3] %f\n",wpl.a[3]);
    printf("a[4] %f\n",wpl.a[4]);
    printf("a[5] %f\n",wpl.a[5]);
    printf("a[6] %f\n",wpl.a[6]);
    printf("a[7] %f\n",wpl.a[7]);
    printf("a[8] %f\n",wpl.a[8]);
    printf("a[9] %f\n",wpl.a[9]);
    printf("a[10] %f\n",wpl.a[10]);

    printf("HOD: %e %e %e %e %e\n",HOD.M_min, HOD.M1, HOD.alpha, alex_mmin, alex_mmax);
    x=log(1.9*pow(3*HOD.M_max/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0));
    printf("%e %e\n",r,x);

    m     = mlo;                // mass in non log units

    r_vir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1
    x     = alex_sat_r/(2*r_vir);  // correct input r/(2*Rvir)
    conc  = halo_concentration(m); 
    fsat  = dFdx_ss(x, conc);      // NFW convolved with itself (EQ 7)
    fcen  = dFdx_cs(x, conc);      // NFW convolved with itself (EQ 7)

    // In thresholds 
    //muh(0);
    //fmuh(m);
    n_sat = Nsat_xigm(m);          // See equation 6 and 7
    n_cen = Ncen_xigm(m);          // See equation 6 and 7

    // check h^1 masses here
    n_halo   = dndM_interp(m);    

    printf("\n mlo %f \n",log10(m));
    printf("r_vir %e \n",r_vir);  
    printf("x %e \n",x);
    printf("conc %e \n",conc);
    printf("fsat %e \n",fsat);
    printf("fcen %e \n",fcen);
    printf("fsat %e \n",fsat);
    printf("n_sat %e \n",n_sat);
    printf("n_cen %e \n",n_cen);

    m     = HOD.M_max;                // mass in non log units

    r_vir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1
    x     = alex_sat_r/(2*r_vir);  // correct input r/(2*Rvir)
    conc  = halo_concentration(m); 
    fsat  = dFdx_ss(x, conc);      // NFW convolved with itself (EQ 7)
    fcen  = dFdx_cs(x, conc);      // NFW convolved with itself (EQ 7)

    // In thresholds 
    //muh(1);
    //fmuh(m);
    n_sat = Nsat_xigm(m);          // See equation 6 and 7
    n_cen = Ncen_xigm(m);          // See equation 6 and 7

    // check h^1 masses here
    n_halo   = dndM_interp(m);    

    printf("\n M_max %f \n",log10(m));
    printf("r_vir %e \n",r_vir);  
    printf("x %e \n",x);
    printf("conc %e \n",conc);
    printf("fsat %e \n",fsat);
    printf("fcen %e \n",fcen);
    printf("fsat %e \n",fsat);
    printf("n_sat %e \n",n_sat);
    printf("n_cen %e \n",n_cen);
    printf("n_halo %e \n",n_halo);
    exit(0);
  }
  return(res);
}

//-----------------------------------------------------------
// Calculate xi for the One Halo For Simple NFW
// This is used for test purposes
//-----------------------------------------------------------

// Input R is in co H=100
// Output xi is in co H=100 
// Output is ALREADY 1+xi
double xi_1h_simple(double r)
{
  double res,norm,ng,mlo,mmin;
  double temp, m,r_vir,x,fsat,fcen,conc,n_sat,n_cen,n_halo,rho_mean;

  m     = 1.461809e12*HUBBLE;                                 // mass in linear units, h inverse
  r_vir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD);  // co h^-1, delta halo from batch file
  x     = r/(2*r_vir);                                        // correct input r/(2*Rvir)
  conc  = halo_concentration(m); 
  fcen  = dFdx_cs(x, conc);                                   // related to NFW 
     
  rho_mean = RHO_CRIT*OMEGA_M;
  
  // Need all the normalization factors (to go from fcen to NFW):
  res = (1.0/(4*PI*pow(r,2)))*(m/(rho_mean))*(1.0/(2*r_vir))*fcen; 

  return(res);
}

//-----------------------------------------------------------
// Function for One Halo (see Eq 6 Yoo et al.)
// Inputs mass in natural log units
//-----------------------------------------------------------

// Add the central and satellite term here together
double funct_1h(double m)
{
  double res, temp, r_vir,x,fsat,fcen,conc,n_sat,n_cen,n_halo,rho_mean;
  int test;
  test=0;

  m     = exp(m);                                            // mass in linear units
  r_vir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1, delta halo from batch file
  x     = alex_sat_r/(2*r_vir);                              // correct input r/(2*Rvir)
  conc  = halo_concentration(m); 
  fsat  = dFdx_ss(x, conc);                                  // related to NFW convolved with itself (EQ 7)
  fcen  = dFdx_cs(x, conc);                                  // related to NFW 

  // In thresholds 
  //muh(2);
  //fmuh(m);
  // printf("FUNC_1H: %e\n",m); 

  n_sat = Nsat_xigm(m);                                      // See equation 6 and 7
  n_cen = Ncen_xigm(m);                                      // See equation 6 and 7
  
  n_halo   = dndM_interp(m);    
  rho_mean = RHO_CRIT*OMEGA_M;

  //  dM * (dn/dM) *(M/rho_mean) * (1/2Rvir) * (n_sat*F_sat + n_cen*F_cent)
  // Add the central and satellite term here

  // Cent and Sat together
  if(ALEX_PRETTY_PLOT_FLAG == 0){
    res = m*n_halo*(m/rho_mean)*(1.0/(2*r_vir))*((n_sat*fsat)+(n_cen*fcen)); // extra m for log integration
  }

  // Only cent
  if(ALEX_PRETTY_PLOT_FLAG == 1){
    res = m*n_halo*(m/rho_mean)*(1.0/(2*r_vir))*(n_cen*fcen); // extra m for log integration
  }

  // Only sat
  if(ALEX_PRETTY_PLOT_FLAG == 2){
    res = m*n_halo*(m/rho_mean)*(1.0/(2*r_vir))*(n_sat*fsat); // extra m for log integration
  }
  
  return(res);
}

//-----------------------------------------------------------
//  ***********   CALCULATE DS   ****************************
//  *********************************************************
//-----------------------------------------------------------

//----------------------------------------------------------- 
// Integrate to get Delta_Sigma
// Input r is in Physical Kpc H=72
//-----------------------------------------------------------

// Input r is in Physical Kpc H=72
// Output Delta Sigma is in physical 0.72
double delta_sigma(double r)
{
  double sigma_bar1, sigma_bar2, r_comov, r_min,x, t1,t2;

  // R is in Co-moving Mpc H=100
  r_comov = r*(1+REDSHIFT)*1e-3*HUBBLE; 

  // Minimum value at which this is evaluated
  r_min=1e-3;

  // Call tabulate_xi
  // In order to avoid nested integrals
  // MUST CALL HERE FIRST
  x = tabulate_xi(0.1);

  // Second Term in Delta Sigma
  sigma_bar2 = xi2sigmabar(r_comov);

  // Mean within sphere
  if(r_comov < r_min){
    printf("exited in delta_sigma: r too small\n");
    printf("r in co-mov Mpc h=1 %f \n",r_comov);
    exit(0);
    }

  if(r_comov >= r_min){
    alex_ds_r = r_comov;
    // !! Be careful of nested integrations here
    // There are static variable in qromo/ call once
    // with mipnt and once with trapzd
    // !! 1e-4 is the minimum for the integration
    sigma_bar1 = qromo(funct1,1e-4,r_comov,alextrapzd);
  }

  return(sigma_bar1-sigma_bar2);
}

// Input r is in Physical Kpc H=72
// Output Sigma is in physical 0.72
double sigma(double r)
{
  double sigma_bar2, r_comov, x;

  // R is in Co-moving Mpc H=100
  r_comov = r*(1+REDSHIFT)*1e-3*HUBBLE; 

  // Call tabulate_xi
  // In order to avoid nested integrals
  // MUST CALL HERE FIRST
  x = tabulate_xi(0.1);

  // Second Term in Delta Sigma
  sigma_bar2 = xi2sigmabar(r_comov);

  return(sigma_bar2);
}

//----------------------------------------------------------- 
// Function for Delta Sigma (calculate mean within sphere
// of radius r)
//-----------------------------------------------------------

// Input r is comoving H=100
// Output Delta Sigma is in physical 0.72
double funct1(double r)
{
  double res;

  res=2.0*r*xi2sigmabar(r)/(alex_ds_r*alex_ds_r);

  return(res);
}

//-----------------------------------------------------------
// From xigm to Sigma_bar
// Integrtion is done in comoving units here
// Need the factor of 2 here!
//-----------------------------------------------------------

// Input r is in Co-moving Mpc H=100
// Output sigmabar is in physical H=72 (conversion in sigmabar_funct)
double xi2sigmabar(double r) 
{
  double res, r_low, r_min;

  r_min=1e-3; 

  if(r < r_min){
    alex_trans_R = r_min;    // Global variable transverse dis
    res = 2*qromo(sigmabar_funct,0,50,alexmidpnt);
  }
  if(r >= r_min){
    alex_trans_R = r;        // Global variable transverse dis
    res = 2*qromo(sigmabar_funct,0,50,alexmidpnt);
   
    if(ERROR_FLAG == 1){
      printf("\n STOPPED IN xi2sigmabar \n",res);
      printf("\n res %f \n",res);
      printf("r %f\n",r);
      exit(0);
    }
  }
  return(res);
}

//-----------------------------------------------------------
// Sigma_Bar Funct
// Function to be integrated to go from xi to Sigma_bar
// input r is in Co-moving Mpc H=100
//-----------------------------------------------------------

// Input r is in Co-moving Mpc H=100
// Output sigmabar is in physical H=72
double sigmabar_funct(double r)
{
  double res, x,r3d;

  // Convert to r3D here
  // R is already in Co-moving Mpc H=100
  // r= r_p and alex_trans_R=r_parallel (alex_trans_R in co h inverse)
  // So the integral is done in comoving 
  r3d = pow(pow(r,2)+pow(alex_trans_R,2),0.5); 

  // I think this is 1+xi not xi .. need to check ...
  x = tabulate_xi(r3d);  // Use tabulated value of xi here

  // To get Sigma in co h inverse only need to multiply by mean density
  //res = x*RHO_CRIT*OMEGA_M; // This is in comov h inverse
  res = x*RHO_CRIT*OMEGA_M; // This is in comov h inverse

  // Now get Sigma in physical h72 [ units are h72 Msun/pc^2]
  res = res*pow((1+REDSHIFT),2)*HUBBLE*1e-12; // this is in physical h72 Msun/pc^2

  //  fprintf(stdout,"%e %e %e %e\n",r,r3d,x,res);
  //fflush(stdout);
  return(res);
}

//-----------------------------------------------------------
//  ***********   TWO HALO  *********************************
//  Includes both Sat and Cent
//  *********************************************************
//-----------------------------------------------------------

// input: h72 physical
// output: h72 physical
void two_halo_ds(double r[8], double *ds)
{
  int i;
  
  ALEX_RESET_FLAG_XI = 1;
  ALEX_XI_OPTION     = 2;     // For 1h 2h (see calc_xi)

  for(i=0;i<=N_LENSING-1;i++){
    ds[i] = delta_sigma(r[i]);
  }
}

//-----------------------------------------------------------
//  ***********   ONE HALO DS   *****************************
//  Includes both Sat and Cent
//  *********************************************************
//-----------------------------------------------------------

// input: h72 physical
// output: h72 physical
void one_halo_ds(double r[8], double *ds)
{
  int i; 

  ALEX_RESET_FLAG_XI = 1;
  ALEX_XI_OPTION     = 1;     // For 1h (see calc_xi)

  for(i=0;i<=N_LENSING-1;i++){
    ds[i] = delta_sigma(r[i]);
  }
}

//-----------------------------------------------------------
//  ***********   FULL MODEL   *****************************
//  Integrates both 2h and 1h together
//  *********************************************************
//-----------------------------------------------------------

// input: h72 physical
// output: h72 physical
void full_mod_ds(double r[8], double *ds)
{
  int i; 

  ALEX_RESET_FLAG_XI = 1;
  ALEX_XI_OPTION     = 3;     // For full model (see calc_xi)

  for(i=0;i<=N_LENSING-1;i++){
    ds[i] = delta_sigma(r[i]);
  }
}

//-----------------------------------------------------------
//  ***********    POINT SOURCE   ***************************
//  *********************************************************
//-----------------------------------------------------------

//-----------------------------------------------------------
//  POINT SOURCE TERM
//  M0    [ h72^-1 Msun  ]
//  X^2   [ h72^-2 Mpc^2 ]   -> h72^1 Msun pc^-2
//-----------------------------------------------------------

// input: h72 physical
// output: h72 physical
void ps_ds(double r[8], double ms, double *ds)
{
  int i;

  for(i=0;i<=N_LENSING-1;i++){
    //Note: convert kpc to Mpc
    ds[i] = (ms/1e12)/(PI*pow(r[i]*1e-3,2));   //factor of 1e12 to convert to pc^2
  }
}

//-----------------------------------------------------------
//  ***********  FULL DS MODEL   ***************************
//  *********************************************************
//-----------------------------------------------------------

// input: h72 physical
// output: h72 physical
void ds_mod(double r[8], double ms, double *ds, double *res_point_source,  double *res_2h, double *res_1h, double *res_full, double *a)
{
  int i,test;
  double t1,t2;
  test=0;

  // ** r is in physical kpc H=72 **

  if(LENSING_OUTPUT_FLAG){
    ps_ds(r,ms,res_point_source);           // Point source
    one_halo_ds(r, res_1h);                 // One halo (including satellites and scatter)
    two_halo_ds(r, res_2h);                 // Two Halo
    full_mod_ds(r, res_full);               // Full model

    for(i=0;i<=N_LENSING-1;i++)
      {
	// Sum of componants: point source + 1h + 2h
	ds[i] = res_point_source[i]+res_1h[i]+res_2h[i];
      
	// Full model should also include point source
	res_full[i] = res_point_source[i]+res_full[i];

	if(test){
	  printf("\n r : %f\n", r[i]);
	  printf("ps : %f\n", res_point_source[i]);
	  printf("1h : %f\n", res_1h[i]);
	  printf("2h : %f\n", res_2h[i]);
	}
      }
  }
  // Integrate everything together
  else{
    ps_ds(r,ms,res_point_source);           // Point source
    //t1=second();
    full_mod_ds(r, res_full);               // Full model, integration using qromo (0.7s*9=6 sec)
    //t2=second();
    //t1=timediff(t1,t2);
    //printf("\n full mod ds_mod : Timing sec: %f\n", t1);

    for(i=0;i<=N_LENSING-1;i++)
      {
	// Sum of componants: point source + full model
	ds[i] = res_point_source[i]+res_full[i];
      
	if(test){
	  printf("\n r : %f\n", r[i]);
	  printf("ps : %f\n", res_point_source[i]);
	  printf("full mod : %f\n", res_full[i]);
	}
      }
  }

  if(test){
    printf("\n EXITING IN CHI2_LENSING,5 \n");
    exit(0);
  }
}

//-----------------------------------------------------------
//  ***********  MISC   *************************************
//  *********************************************************
//-----------------------------------------------------------

//-----------------------------------------------------------
//  STELLAR MASS TO HALO MASS
//  ms is stellar mass in non log units
//-----------------------------------------------------------

// input: h72 
// output: h72 
double ms_to_mhalo(double ms, double *a)
{
  double res, mhalo_norm, mstellar_norm,beta,delta,gamma,x;
  int test;
  test=0;
  res=0.0;

  mhalo_norm     = a[1];          // M1 (in LOG10 units)   
  mstellar_norm  = pow(10,a[2]);  // M*,0 (convert to linear units ...)
  beta           = a[3];          // beta
  delta          = a[4];          // delta
  gamma          = a[5];          // gamma

  x=ms/mstellar_norm;             // in lin units

  // This is in log10
  res = mhalo_norm+(beta*log(x)/log(10))+(pow(x,delta)/(1+pow(x,-gamma)))-0.5;

  // Now convert to linear units
  res = pow(10,res);

  if(test){
    printf("\n");
    printf("This is the stellar mass : %e\n", log(ms)/log(10));
    printf("%f %f %f %f %f\n",a[1],a[2],a[3],a[4],a[5]);
    printf("This is the halo mass given the stellar mass : %e\n", log(res)/log(10));
  }
  return(res);
}

//-----------------------------------------------------------
// Make a pretty plot for the paper
//-----------------------------------------------------------
 
void make_pretty_plot(double ms, int plotnum){
 
  double r[1000];
  char fname[1000], aa[1000];
  FILE *fp;
  int i;
  double ps[1000],all_1h[1000],sat_1h[1000],all_2h[1000],sm_lim_0,sm_lim_1,mass[500],n_sat,n_cen;


  // make a file with r
    for(i=0;i<=999;i++){
    r[i]  = (i*8+10);   // in kpc h72

    // input: h72 physical kpc
    ps[i] = (ms/1e12)/(PI*pow(r[i]*1e-3,2));   //factor of 1e12 to convert to pc^2
  }

  // 1h cent
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=0;
  for(i=0;i<=999;i++){
    all_1h[i] = delta_sigma(r[i]);
  }  
  
  // 1h sat
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=2;
  for(i=0;i<=999;i++){
    sat_1h[i] = delta_sigma(r[i]);
  } 
  
  // 2h
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=2;
  for(i=0;i<=999;i++){
    all_2h[i] = delta_sigma(r[i]);
  } 

  if(plotnum == 1){
    fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_1.txt","w");
  }
  if(plotnum == 2){
    fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_2.txt","w");
  }

  // Write the output here
  for(i=0;i<=999;i++){
    fprintf(fp,"%lf %lf %lf %f %f\n",r[i],ps[i],all_1h[i],sat_1h[i],all_2h[i]);
  }
  fclose(fp);

  // NOW MAKE A FILE THAT CONTAINS THE HOD
  // Mass in non log units
  for(i=0;i<=500;i++){
    mass[i]  = pow(10.0, 10+(i*0.01));
  }

  /*this_sm1  = 10.0^9.0   ;1
    this_sm2  = 10.0^9.5   ;2
    this_sm3  = 10.0^10.0  ;3
    this_sm4  = 10.0^10.5  ;4
    this_sm5  = 10.0^11.0  ;5
    this_sm6  = 10.0^11.5  ;6*/

  //---------
  ALEX_RESET_FLAG_XI=1;
  sm_lim_0 = pow(10,  9.0);
  sm_lim_1 = pow(10,  9.5);         
  set_up_hod_for_lensing(sm_lim_0,sm_lim_1, wpl.a);      // This now sets up the HOD in thresholds correctly  
  fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_hod_model_1.txt","w");
  for(i=0;i<=500;i++){
    n_sat = Nsat_xigm(mass[i]);          // See equation 6 and 7
    n_cen = Ncen_xigm(mass[i]); 
    fprintf(fp,"%lf %lf %lf \n",log10(mass[i]),n_sat,n_cen);
  }
  fclose(fp);

  //---------
  ALEX_RESET_FLAG_XI=1;
  sm_lim_0 = pow(10,  9.5);
  sm_lim_1 = pow(10,  10.0);
  set_up_hod_for_lensing(sm_lim_0,sm_lim_1 , wpl.a);      // This now sets up the HOD in thresholds correctly  
  fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_hod_model_2.txt","w");
  for(i=0;i<=500;i++){
    n_sat = Nsat_xigm(mass[i]);          // See equation 6 and 7
    n_cen = Ncen_xigm(mass[i]); 
    fprintf(fp,"%lf %lf %lf \n",log10(mass[i]),n_sat,n_cen);
  }
  fclose(fp);

  //---------
  ALEX_RESET_FLAG_XI=1;
  sm_lim_0 = pow(10,  10.0);
  sm_lim_1 = pow(10,  10.5);
  set_up_hod_for_lensing(sm_lim_0,sm_lim_1 , wpl.a);      // This now sets up the HOD in thresholds correctly  
  fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_hod_model_3.txt","w");
  for(i=0;i<=500;i++){
    n_sat = Nsat_xigm(mass[i]);          // See equation 6 and 7
    n_cen = Ncen_xigm(mass[i]); 
    fprintf(fp,"%lf %lf %lf \n",log10(mass[i]),n_sat,n_cen);
  }
  fclose(fp);

  //---------
  ALEX_RESET_FLAG_XI=1;
  sm_lim_0 = pow(10,  10.5);
  sm_lim_1 = pow(10,  11.0);
  set_up_hod_for_lensing(sm_lim_0,sm_lim_1 , wpl.a);      // This now sets up the HOD in thresholds correctly  
  fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_hod_model_4.txt","w");
  for(i=0;i<=500;i++){
    n_sat = Nsat_xigm(mass[i]);          // See equation 6 and 7
    n_cen = Ncen_xigm(mass[i]); 
    fprintf(fp,"%lf %lf %lf \n",log10(mass[i]),n_sat,n_cen);
  }
  fclose(fp);

  //---------
  ALEX_RESET_FLAG_XI=1;
  sm_lim_0 = pow(10,  11.0);
  sm_lim_1 = pow(10,  11.5);
  set_up_hod_for_lensing(sm_lim_0,sm_lim_1 , wpl.a);      // This now sets up the HOD in thresholds correctly  
  fp=fopen("/Users/alexie/Work/SM_paper/pretty_plot_hod_model_5.txt","w");
  for(i=0;i<=500;i++){
    n_sat = Nsat_xigm(mass[i]);          // See equation 6 and 7
    n_cen = Ncen_xigm(mass[i]); 
    fprintf(fp,"%lf %lf %lf \n",log10(mass[i]),n_sat,n_cen);
  }
  fclose(fp);

  printf("End of pretty plot");

}

//-----------------------------------------------------------
// Sigma for Fabian
//-----------------------------------------------------------
 
void sigma_fabian(double ms){
 
  double r[1000];
  char fname[1000], aa[1000];
  FILE *fp;
  int i;
  double ps[1000],all_1h[1000],sat_1h[1000],all_2h[1000],full_mod[1000],sm_lim_0,sm_lim_1;

  // make a file with r
  for(i=0;i<=999;i++){
    r[i]  = (i*8+10);   // in kpc h72
  }

  // *** DONT USE TABULATION OF 2H HERE:
  ALEX_USE_2H_TABULATE=0;

  // 1h cent
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=0;
  for(i=0;i<=999;i++){
    all_1h[i] = sigma(r[i]);
  }  
  
  // 1h sat
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=2;
  for(i=0;i<=999;i++){
    sat_1h[i] = sigma(r[i]);
  } 
  
  // 2h
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=2;
  for(i=0;i<=999;i++){
    all_2h[i] = sigma(r[i]);
  } 

  // full mode
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=3;
  for(i=0;i<=999;i++){
    full_mod[i] = sigma(r[i]);
  } 

  fp=fopen("/Users/alexie/Work/Magnification/Fabian/for_eric_11.0_12.0.txt","w");
  //fp=fopen("/Users/alexie/Work/Magnification/Fabian/for_eric_10.8_12.0.txt","w");

  // Write the output here
  for(i=0;i<=999;i++){
    fprintf(fp,"%lf %lf %lf %lf %lf \n",r[i],all_1h[i],sat_1h[i],all_2h[i],full_mod[i]);
  }
  fclose(fp);

  exit(0);
  printf("End of sigma_fabian");
}

//-----------------------------------------------------------
// 2h term for EDO
//-----------------------------------------------------------
 
void make_pretty_plot_foredo(double ms){
 
  double r[1000];
  char fname[1000], aa[1000];
  FILE *fp;
  int i;
  double ps[1000],all_1h[1000],sat_1h[1000],all_2h[1000], m1;

  // make a file with r
  for(i=0;i<=999;i++){
    r[i]  = (i*8+10);   // in kpc h72
  }
  
  // 2h
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=2;
  for(i=0;i<=999;i++){
    all_2h[i] = delta_sigma(r[i]);
  } 

  fp=fopen("/Users/alexie/Work/SM_paper/alexie_2h_gg.txt","w");
  
  // Write the output here
  for(i=0;i<=999;i++){
    fprintf(fp,"%lf %lf\n",r[i],all_2h[i]);
  }
  printf("End of pretty plot");
  fclose(fp);
}

// This make a file of mass, concentration, and z that
// I use to make an IDL version of the mass-concentraion
// relation for mass definition = 200 * back

void write_text_halo_concentration(void){
 
  double log_mass[400];
  char fname[400], aa[400];
  FILE *fp;
  int i,j;
  double conc;

  // make a file with mass
  // from 7 to 16.5
  // in log10 units
  for(i=0;i<=399;i++){
    log_mass[i]  = (i*((16.5-7.0)/399.0)+7); 
  }
  
  fp=fopen("/Users/alexie/idl/Conc_munoz/conc_munoz.txt","w");

  REDSHIFT=0.0;
  RESET_COSMOLOGY++;

  // Write the output here
  for(j=0;j<=10;j++){
    for(i=0;i<=399;i++){
      conc = halo_concentration(pow(10.0,log_mass[i]));
      fprintf(fp,"%lf %lf %lf\n",log_mass[i],REDSHIFT, conc);
    }
    REDSHIFT=REDSHIFT+0.1;
    RESET_COSMOLOGY++;
  }

  fclose(fp);
  printf("\n End of write_text_halo_concentration \n");
  exit(0);
}

//-----------------------------------------------------------
//  ***********    TEST FUNCTIONS   *************************
//  *********************************************************
//-----------------------------------------------------------

void test_delta_sigma(void){
 
  // Need to set scatter to small values to use this code :
  double testr[100];
  char fname[100], aa[100];
  FILE *fp;
  int i,j;
  double temp,sm0,mh,temp2,conc,r_vir,x,fcen,r_vir_phy,rho_nfw;

  //sm0 = pow(10.0,11.2); -> 1.090702e13
  //printf("ms to mhalo : %e \n",ms_to_mhalo(sm0,wpl.a));
  //exit(0);

  //----------------------------------
  // The following is used in plot_covariance2.pro

  RESET_COSMOLOGY++;

  ALEX_XI_TEST_OPTION=0; // Uses xi_1h_simple to test NFW
  ALEX_2H_BIN_FLAG=0;
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=1; // 1h cent (NFW ONLY)

  REDSHIFT=0.372; 
  alex_mmin = pow(10.0,11.225);           // SM bins
  alex_mmax = pow(10.0,11.229);           // SM bins
  sm0 = pow(10.0,11.228);

  printf("\n \n redshift %f \n :",REDSHIFT);
  mh=ms_to_mhalo(sm0,wpl.a); // h72
  printf("ms to mhalo (h72 units) : %e \n",mh);
  printf("ms to mhalo (h inverse units) : %e \n",mh*HUBBLE);
  conc=halo_concentration(mh*HUBBLE);
  printf("conc %f \n", conc);

  r_vir = pow(3*mh*HUBBLE/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1, delta halo from batch file
  r_vir_phy =   r_vir*(1.0/(1+REDSHIFT))*1.0/HUBBLE;
  printf("rvir comov h^-1: %e \n",r_vir);
  printf("rvir phy h72: %e \n",r_vir_phy);
  printf("DELTA HALO : %f \n",DELTA_HALO);

  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly

  // make a file with r, in phy h72
  // This should be in kpc
  for(i=0;i<=99;i++){
    testr[i]  = (i*((1000.0-20)/99.0)+20); 
  }
  
  fp=fopen("/Users/alexie/Work/HOD/Pro/z0.372_m11.22.txt","w");

  // NOTE: rho_nfw (r3d) = (1/4pi*r^2)*M*(1/2*rvir)*F'(x)
  // Write the output here
  for(i=0;i<=99;i++){
    temp = delta_sigma(testr[i]); 
    temp2 = sigma(testr[i]); 
    x     = testr[i]*1e-3/(2*r_vir_phy);
    fcen  = dFdx_cs(x, conc);
    rho_nfw = (1.0/(4*PI*pow(testr[i]*1e-3,2)))*mh*(1.0/(2*r_vir_phy))*fcen; // in phy h72
    fprintf(fp,"%lf %lf %lf %lf\n",testr[i], temp, temp2,rho_nfw);
  }
     
  fclose(fp);

  //----------------------------------
  // The following are used for test_chi2_lensing.pro

  RESET_COSMOLOGY++;

  ALEX_XI_TEST_OPTION=0; // Uses xi_1h_simple to test NFW
  ALEX_2H_BIN_FLAG=0;
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=1; // 1h cent (NFW ONLY)

  REDSHIFT=0.3; 
  alex_mmin = pow(10.0,10.795);           // SM bins
  alex_mmax = pow(10.0,10.805);           // SM bins
  sm0 = pow(10.0,10.8);

  printf("\n \n redshift %f \n :",REDSHIFT);
  mh=ms_to_mhalo(sm0,wpl.a); // h72
  printf("ms to mhalo (h72 units) : %e \n",mh);
  printf("ms to mhalo (h inverse units) : %e \n",mh*HUBBLE);
  conc=halo_concentration(mh*HUBBLE);
  printf("conc %f \n", conc);

  r_vir = pow(3*mh*HUBBLE/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1, delta halo from batch file
  r_vir_phy =   r_vir*(1.0/(1+REDSHIFT))*1.0/HUBBLE;
  printf("rvir comov h^-1: %e \n",r_vir);
  printf("rvir phy h72: %e \n",r_vir_phy);
  printf("DELTA HALO : %f \n",DELTA_HALO);

  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly

  // make a file with r, in phy h72
  // This should be in kpc
  for(i=0;i<=99;i++){
    testr[i]  = (i*((1000.0-20)/99.0)+20); 
  }
  
  fp=fopen("/Users/alexie/Work/HOD/Pro/z0.3_m10.8.txt","w");

  // NOTE: rho_nfw (r3d) = (1/4pi*r^2)*M*(1/2*rvir)*F'(x)
  // Write the output here
  for(i=0;i<=99;i++){
    temp = delta_sigma(testr[i]); 
    temp2 = sigma(testr[i]); 
    x     = testr[i]*1e-3/(2*r_vir_phy);
    fcen  = dFdx_cs(x, conc);
    rho_nfw = (1.0/(4*PI*pow(testr[i]*1e-3,2)))*mh*(1.0/(2*r_vir_phy))*fcen; // in phy h72
    fprintf(fp,"%lf %lf %lf %lf\n",testr[i], temp, temp2,rho_nfw);
  }
     
  fclose(fp);

  //----------------------------------

  RESET_COSMOLOGY++;

  ALEX_XI_TEST_OPTION=0; // Uses xi_1h_simple to test NFW
  ALEX_2H_BIN_FLAG=0;
  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1; 
  ALEX_PRETTY_PLOT_FLAG=1; // 1h cent (NFW ONLY)

  REDSHIFT=0.7; 
  alex_mmin = pow(10.0,10.795);           // SM bins
  alex_mmax = pow(10.0,10.805);           // SM bins
  sm0 = pow(10.0,10.8);

  printf("\n \n redshift %f \n :",REDSHIFT);
  mh=ms_to_mhalo(sm0,wpl.a); // h72
  printf("ms to mhalo (h72 units) : %e \n",mh);
  printf("ms to mhalo (h inverse units) : %e \n",mh*HUBBLE);
  conc=halo_concentration(mh*HUBBLE);
  printf("conc %f \n", conc);

  r_vir = pow(3*mh*HUBBLE/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD); // co h^-1, delta halo from batch file
  r_vir_phy =   r_vir*(1.0/(1+REDSHIFT))*1.0/HUBBLE;
  printf("rvir comov h^-1: %e \n",r_vir);
  printf("rvir phy h72: %e \n",r_vir_phy);
  printf("DELTA HALO : %f \n",DELTA_HALO);

  set_up_hod_for_lensing(alex_mmin, alex_mmax, wpl.a);      // This now sets up the HOD in thresholds correctly

  // make a file with r, in phy h72
  // This should be in kpc
  for(i=0;i<=99;i++){
    testr[i]  = (i*((1000.0-20)/99.0)+20); 
  }
  
  fp=fopen("/Users/alexie/Work/HOD/Pro/z0.7_m10.8.txt","w");

  // NOTE: rho_nfw (r3d) = (1/4pi*r^2)*M*(1/2*rvir)*F'(x)
  // Write the output here
  for(i=0;i<=99;i++){
    temp = delta_sigma(testr[i]); 
    temp2 = sigma(testr[i]); 
    x     = testr[i]*1e-3/(2*r_vir_phy);
    fcen  = dFdx_cs(x, conc);
    rho_nfw = (1.0/(4*PI*pow(testr[i]*1e-3,2)))*mh*(1.0/(2*r_vir_phy))*fcen; // in phy h72
    fprintf(fp,"%lf %lf %lf %lf\n",testr[i], temp, temp2,rho_nfw);
  }
     
  fclose(fp);

  //----------------------------------

  printf("\n >> End of test_delta_sigma \n");
  exit(0);
}

// Test xi, sigma_bar, and delta_sigma
void test_xi(double n)
{
  double r[1000],x[1000],s[1000],ds[1000],fsat[1000],fcen[1000],r_low;
  int i;
  char fname[1000], aa[1000];
  FILE *fp;
  double m,r_vir,xx,conc;

  printf("\n --------- TEST XI ---------------\n");

  ALEX_RESET_FLAG_XI=1;
  ALEX_XI_OPTION=1;            
  fp=fopen("/Users/alexie/Work/Jeremy/Test/test_xi.txt","w");

  for(i=0;i<=999;i++){
    r[i]  = (i*5+1)*1e-3;   // comoving Mpc
    s[i]  = xi2sigmabar(r[i]);
    //printf("\n >>>>> ALEX_RESET_FLAG %d:\n",ALEX_RESET_FLAG_XI);
    x[i]  = tabulate_xi(r[i]);
    //sigma_bar2 = xi2sigmabar(r_comov);
    // Note: ds must be called in real kpc
    //ds[i] = delta_sigma(r[i]/((1+REDSHIFT)*1e-3*HUBBLE)); 
    ds[i] = delta_sigma(r[i]);
    fprintf(fp,"%lf %lf %lf %f\n",r[i],x[i],s[i],ds[i]);}

  fclose(fp);
 
  printf("\n >>>>> \n");

  printf("\n --------- EXIT TEST XI ---------------\n");
  exit(0);
}

// THIS PROGRAM IS NOT FINISHED YET ...
// See if the code can correctly predict the PDF of Sigma
void write_sigma_pdf(double r[8]){
 
  double log_mass[400];
  char fname[400], aa[400];
  FILE *fp;
  int i,j;
  double conc;

  // make a file with mass
  // in log10 units
  for(i=0;i<=399;i++){
    log_mass[i]  = (i*((16.0-9.0)/399.0)+9); 
  }
  
  fp=fopen("/Users/alexie/Work/HOD/Distributions/sm0_r0_distribution.txt","w");

  //****

  // Integrate Fcen(x) and Fcen(x) 
  // Use this to predict Sigma(m,r)
  // then check this result
  

  // Write the output here
  for(j=0;j<=10;j++){
    for(i=0;i<=399;i++){
      conc = halo_concentration(pow(10.0,log_mass[i]));
      fprintf(fp,"%lf %lf %lf\n",log_mass[i],REDSHIFT, conc);
    }
    REDSHIFT=REDSHIFT+0.1;
    RESET_COSMOLOGY++;
  }

  fclose(fp);
  printf("\n End of write_text_halo_concentration \n");
  exit(0);
}
