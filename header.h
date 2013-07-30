#include "nrutil.h"
#include "stdlib.h"
#include "stdio.h"


/* Function prototypes--> SHMR
 */
double chi2_lensing(double *a);
double every_fucking_observers_mass_function(double mass);
double xi_2h_gm(double r, double mlo, double mhi, double *a);
double Nsat_xigm(double m);
double Ncen_xigm(double m);
double func_galaxy_density_xigm(double m);
double integrate(double (*fxn)(double), double a, double b,int closed_a,int closed_b);
double fsigma(double sig);
double bias_fsigma(double sig, double r);
double sigma2mass(double sigma);
double Ncen_CLF(double m);
double ms_to_mhalo_inversion(double mh);
double ms_to_mhalo(double mh, double *a);
void shmr_clustering(void);
void shmr_color_clustering(void);
void input_stellar_mass_bins(void);
int set_up_hod_for_shmr(double mlo, double mhi, double *a);
void shmr_lensing(void);
void shmr_color_lensing(void);
double pca_chi2(double *xm, int ii);
void chain_postprocess(void);

/* Function prototypes--> utility files
 */
double second(void);
double timediff(double t0,double t1);
void endrun(char *instring);
FILE *openfile(char *ff);
int filesize(FILE *fp);
void least_squares(double *x, double *y, int n, double *a, double *b);
void check_for_smoothness(double *x, double *y, int n, double r);
void read_parameter_file(char *fname);
void output_parameter_file(char *fname);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

/* Function prototypes--> general control of tasks
 */
void tasks(int argc, char **argv);

/* Function prototypes--> numerical recipes
 */
void sort(unsigned long n, double arr[]);
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double midpnt(double (*func)(double), double a, double b, int n);
double alexmidpnt(double (*func)(double), double a, double b, int n);
double alex1midpnt(double (*func)(double), double a, double b, int n);
double alex2midpnt(double (*func)(double), double a, double b, int n);
double alex3midpnt(double (*func)(double), double a, double b, int n);
double midinf(double (*funk)(double), double aa, double bb, int n);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
double zbrent(double (*func)(double), double x1,double x2, double tol);
double qtrap(double (*func)(double), double a, double b, double EPS);
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	    double (*func)(double []));
void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk);
void gaussj(double **a, int n, double **b, int m);
float gasdev(long *idum);
void jacobi(double **a, int n, double d[], double **v, int *nrot);
float ran1(long *idum);
double ran2(long *idum);
double trapzd(double (*func)(double), double a, double b, int n);
double alextrapzd(double (*func)(double), double a, double b, int n);
void sort2(unsigned long n, float arr[], int id[]);
void sort2dbl(unsigned long n, double arr[], int id[]);


/* Function prototypes--> power spectrum routines.
 */
double transfnc(double xk);
double sigmac(double rad);
double nonlinear_sigmac(double rad);
double transfunc_file(double xk);
double nonlinear_power_spectrum(double xk);
double linear_power_spectrum(double xk);
double mstar(void);
double tf_eisenstein_hu(double k);
double cobenorm(double omega_m);
double cobe_prior(double omega_m);
double sigmac_interp(double m);
double sigmac_radius_interp(double m);

/* Function prototypes--> matter correlation function routines.
 */
double xi_interp(double r);
double xi_linear_interp(double r);
double projected_xi_matter(double r);

/* Function protoypes--> outputting matter stats
 */
void output_matter_power_spectrum(void);
void output_matter_correlation_function(void);
void output_matter_variance(void);
void output_halo_concentrations(void);
void output_halo_mass_function(void);
void output_halo_correlation_function(double mass);

/* Function prototypes--> HOD functions (+related);
 */
double N_avg(double m);
double N_cen(double m);
double N_sat(double m);
double moment_ss(double m);
double func_galaxy_density(double m);
double func_central_density(double m);
double func_central_bias(double m);
double func_satfrac(double m);
void set_HOD_params(void);
double set_low_mass(void);
double set_high_central_mass(void);
double N_cen_i(double m, int ii);
double satellite_blue_fraction(double m);
double central_blue_fraction(double m);
double func_halo_density(double m);

double func_central_density(double m);
double func_satellite_density(double m);

double number_weighted_central_mass(void);
double number_weighted_halo_mass(void);

/* Density-dependent HOD functions
 */
void dd_hod_functions(float *hmass, float *hdensity, int nhalo);
int internal_populate_simulation(float *x, float *y, float *z, float reduction_factor, int ivel, 
				 float *vx, float *vy, float *vz);


/* Function prototypes--> HOD functions for second HOD (+related);
 */
double N_avg2(double m);
double N_cen2(double m);
double N_sat2(double m);
double moment_ss2(double m);
double func_galaxy_density2(double m);
double func_satfrac2(double m);
void set_HOD2_params(void);
double set_low_mass2(void);
double set_high_central_mass2(void);

/* Function prototypes--> halo mass function.
 */
double halo_mass_function(double mass);
double dndM_interp(double m);

/* Function prototypes--> halo bias.
 */
double bias_interp(double m, double r);
double bias(double m);
double func_galaxy_bias(double m);

/* Function prototypes--> halo profile: pair density
 */
double dFdx_ss(double x, double c_NFW);
double dFdx_cs(double x, double c_NFW);
double dFdx_ss_interp(double r, double c);
double dFdx_cs_interp(double r, double c);
double nfw_transform(double xk, double m);

/* Function prototypes--> galaxy-matter cross-correlation function
 */
double one_halo_galaxy_matter(double r);
double two_halo_galaxy_matter(double r);
double delta_sigma(double rr);

/* Function prototypes--> real-space galaxy correlation function.
 */
double one_halo_real_space(double r);
double two_halo_real_space(double r);
double ellipsoidal_exclusion_probability(double rv, double r);
double restricted_number_density(double r);
double projected_xi(double r);
double projected_xi_1halo(double r);
double projected_xi_rspace(double r);
double nbody_xi(double r);
double integrated_wp_bin(double r);
double wtheta(double theta);

/* Function prototypes--> redshift-space galaxy correlation function.
 */
double one_halo(double rs, double rp);
double two_halo(double rs, double rp);
void xi_multipoles(void);
double small_scale_measure(double rs);
void calc_rhalf(double r[], double rhalf[], int nr);
double integrated_bin(double xlo, double ylo, double dx1, double dy1, int n);
void linlin_bins(void);
double xi2d_interp(double rs1, double rp1, double rs2, double rp2);
double xi2d_interp_polar(double rs1, double rs2, double phi1, double phi2);


/* Function prototypes--> halo pairwise velocity model
 */
double spherical_collapse_model(double delta);
double galaxy_prob_vz(double vel, double rad, double theta);
void vdelta_v4(double m0, double m1, double rad, double theta, double binsize, 
	       double v0, int nv, double *a, double *pv, double *pt, double *pz,
	       double wgal[3], double sgal[4]);
void pairwise_velocity_dispersion(void);
void output_velocity_moments(int imass);
double jeans_dispersion(double mass, double rx, double v[]);


/* Function prototypes--> halo concentrations.
 */
double growthfactor(double z);
double halo_concentration(double m);
double cvir_model(double mass);
double munoz_cuartas_cvir_model(double mass);
double halo_mass_conversion(double mvir, double *cvir1, double delta_halo);
double halo_mass_conversion2(double mvir, double cvir1, double delta_vir, double delta_halo);
double halo_c200(double m);
double HK_func(double x);

/* Function prototypes--> HOD fitting of correlation functions.
 */
void shmr_minimization(void);
void wp_minimization(char *fname);
void mcmc_minimization(void);
void zspace_minimization(char *fname);
double chi2_wp(double *a);
double chi2_zspace(double *a);
void fit_color_samples(void);

/* Function prototypes--> Linear+exponential model.
 */
double kaiser_distortion(double rs, double rp);
void fit_dispersion_model(void);
double linear_kaiser_distortion(double rs, double rp);

/* Function prototypes--> Jean's solution to velocity dispersions.
 */
double satellite_dispersion(double r, double rs, double rvir);
double velocity_dispersion_profile(double sig, double cvir, double rvir, double r);

/* Function prototypes--> random position/velocities of sat gals (simulation population)
 */
void populate_simulation(void);
void populate_sampled_simulation(void);
double NFW_density(double r, double rs, double ps);
double NFW_velocity(double mass, double v[], double mag);
double NFW_position(double mass, double x[]);
int poisson_deviate(double nave);
double poisson_prob(int n, double nave);

/* Definitions
 */
#define PI       3.14159265358979323846
#define TWOPI    (2*PI)
#define THIRD    (1.0/3.0)
#define ROOT2    1.41421356237309504880
#define RT2PI    2.50662827463100050241
#define LN_2     0.6931471805599452
#define ROOT8    2.82842712475
#define WORKBUF  1000
#define LOGE_10  2.30258509
#define c_on_H0  3000

#define mabs(A)  ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define muh(x)   fprintf(stdout,"%d\n",x);fflush(stdout)
#define fmuh(x)  fprintf(stdout,"%e\n",x);fflush(stdout)
#define square(x) (x*x)

/* Global variables
 */
extern double 
  GAMMA,
  HUBBLE,
  HUBBLE_UNITS,
  HUBBLEZ,
  SIGMA_8,
  SIGMA_8Z0,
  RHO_CRIT,
  SPECTRAL_INDX,
  OMEGA_M,
  OMEGA_Z,
  OMEGA_B,
  OMEGA_TEMP,
  DELTA_CRIT,
  MSTAR,
  GALAXY_DENSITY,
  GALAXY_DENSITY2,
  NG_ONE_HALO,
  GRAVITY,
  VZ_LIMIT,
  BOX_SIZE,
  RESOLUTION,
  R_MIN_2HALO,
  MASS_PER_PARTICLE,
  GALAXY_BIAS,
  DELTA_HALO,
  VBIAS,
  VBIAS_SLOPE,
  VBIAS_C,
  VBIAS_MASS_THRESHOLD,
  CVIR_FAC,
  MAXVEL,
  BIAS_A,
  BIAS_B,
  BIAS_C,
  JENKINS_A,
  JENKINS_B,
  JENKINS_C,
  DNDM_PARAMS[10],
  SIGV,
  BETA,
  XI_MAX_RADIUS,
  MASS_THRESHOLD,
  DENSITY_THRESHOLD,
  REDSHIFT,
  LOCAL_DENSITY,
  EXCLUSION_RADIUS,
  SATELLITE_FRACTION;

extern int RESET_COSMOLOGY,
  RESET_KAISER,
  ITRANS,
  COVAR,
  PCA,
  COVARZ,
  EXCLUSION,
  SOFT_CENTRAL_CUTOFF,
  NUM_POW2MASS_BINS,
  FIX_PARAM,
  DEPROJECTED,
  OUTPUT,
  MCMC_OUTPUT,
  POWELL,
  MCMC,
  BEST_FIT,
  LINEAR_PSP,
  KAISER,
  RESET_PVZ,
  ERROR_FLAG,
  RESET_FLAG_2H,
  RESET_FLAG_1H,
  GAO_EFFECT,
  IVFLAG,
  WP_ONLY,
  N_HOD_PARAMS,
  XCORR,
  RESTART,
  USE_ERRORS,
  DENSITY_DEPENDENCE;

extern char RESTART_FILE[100];

extern long IDUM_MCMC;

extern int ARGC;
extern char **ARGV;

extern int *used_halo;

/* Stuff for wtheta
 */
extern double MSHIFT,
  MALL,
  NGAL_DRG;

/* These are variables for the parallel implementation.
 */
extern int ThisTask,NTask;


/* HOD parameters. For definitions, look at the comments in
 * hod_functions.c
 */
extern struct hod_parameters {
  double M_min;
  double M_max;
  double M1;
  double M_cut;
  double M_low,M_low0;
  double M_hi;
  double sigma_logM;
  double alpha;
  double MaxCen;
  double M_cen_max;

  double fblue0_cen;
  double sigma_fblue_cen;
  double fblue0_sat;
  double sigma_fblue_sat;
  double blue_fraction;

  int color;
  int pdfc;
  int pdfs;
  int free[100];
  int i_wp;

  double M_sat_break;
  double alpha1;

  double M_min_loden;
  double M_min_hiden;
  double M_min_fac;

  double fredc, freds;
  double mass_shift, mshift2;
  double shift_alpha;

  // parameters for the stellar-to-halo mass relation SHMR
  double shmr[20];
  double red_central_fraction[10];
  double red_central_halomass[10];
  int red_central_nbin;

} HOD, HOD2, HODt;


/* Structure to keep information/data about fitting
 * color-defined samples.
 */
struct COLOR_DATA {
  int ON;
  double ngal_red;
  double ngal_blue;
  double ngal_full;
  int n_red;
  int n_blue;
  int n_full;
  double *r_red,*r_blue,*r_full;
  double *e_red,*e_blue,*e_full;
  double *x_red,*x_blue,*x_full;
  double **covar_red, **covar_blue, **covar_full;
} wp_color;


/* This is to put the work done in xi_multipoles into 
 * a global space.
 */
extern struct workspace {
  int nrad;
  int n_mono;
  int n_quad;
  int n_half;
  int n_z;
  int SDSS_bins;
  double rad[WORKBUF];
  double r_mono[WORKBUF];
  double r_quad[WORKBUF];
  double r_half[WORKBUF];
  double rsigma[WORKBUF][WORKBUF];
  double rpi[WORKBUF][WORKBUF];
  double xi_mono[WORKBUF];
  double xi_quad[WORKBUF];
  double xi_half[WORKBUF];
  double xi_z[WORKBUF][WORKBUF];
  double data_m[WORKBUF];
  double covar_m[WORKBUF][WORKBUF];
  double data_q[WORKBUF];
  double covar_q[WORKBUF][WORKBUF];
  double data_h[WORKBUF];
  double covar_h[WORKBUF][WORKBUF];
  double data_z[WORKBUF][WORKBUF];
  double err_m[WORKBUF];
  double err_q[WORKBUF];
  double err_h[WORKBUF];
  double err_z[WORKBUF][WORKBUF];

  double **covarz;
  int n_ze;
  char covarz_file[WORKBUF];

  double rmlo,rmhi;
  double rqlo,rqhi;
  double rhlo,rhhi;
  double zlo,zhi;

  int ncf;
  double cf[WORKBUF];
  int chi2;
  int imono;
  int imono_covar;
  int iquad;
  int iquad_covar;
  int ihalf;
  int ihalf_covar;
  int izspace;
  int imodel;
  int i_quad2mono;
  int i_monopole;
  char monofile[WORKBUF];
  char monocovarfile[WORKBUF];
  char quadfile[WORKBUF];
  char quadcovarfile[WORKBUF];
  char halffile[WORKBUF];
  char halfcovarfile[WORKBUF];
  char esysfile[WORKBUF];
  char zfile[WORKBUF];
  double absolute_error;
  double percentage_error;
  double eq_abs;
  double em_abs;
  double em_per;
  double eh_per;
  int use_asymptotic_values;
  int iout;
  int SysErrFlag;

  double *psp[6];
  double r_2h[6][100];
  double xi_2h[6][100];
  int n_2h[6];

} Work;

/* Various input files and flags on whether or not to use them.
 */
extern struct file_parameters {
  char HaloFile[1000];
  char HaloDensityFile[1000];
  char TF_file[100];
  char TwoHaloFile[100];
  int  i_TwoHalo;
  char MassFuncFile[100];
  int i_MassFunc;
  char PDFTable[100];
  int i_PDFTable;
  char PSPFile[100];
  int i_Cvir;
  char CvirFile[100];

  // files for SHMR stuff
  int UseStellarMassBinsClustering;
  int UseStellarMassBinsLensing;
  char StellarMassBinsClustering[1000];
  char StellarMassBinsLensing[1000];
  char ClusteringFilenames[1000];
  char ClusteringCovarFilenames[1000];
  char StellarMassFunctionFilename[1000];
  char StellarMassFunctionCovarFilename[1000];
  char LensingFilenames[1000];
  char MCMCRestartFile[1000];
} Files;

/* Various tasks the the program will perform
 */

extern struct perform_tasks {
  int clustering;
  int lensing;
  int wp_minimize;
  int shmr_minimize;
  int MCMC;
  int HOD;
  int populate_sim;
  int matter_xi;
  int matter_pk;
  int sigma_r;
  int cvir;
  int halostats;
  int SHMR;
  int fit_for_bias;
  char root_filename[100];
} Task;


/* Workspace for w_p minimzation.
 */
struct WP {
  double **covar;
  int np;
  int ncf;
  int ncf_tot;
  double pi_max;
  double *r;
  double *x;
  double *e;
  char fname_covar[100];
  char fname_wp[100];
  int format;
  int iter;
  double esys;
  double *eigen;
  int npca;
  double ngal;
  double ngal_err;
  int n_wp;
  
  double fsat_all, fsat_red, fsat_blue;
} wp;

/* Data structure for wp analysis in LENSING with ALEXIE
 */
struct WPL {
  // commands for running the chains
  double stepfac,
    stepfac_burn; 
  int dont_fit_lensing,
    dont_fit_clustering,
    dont_fit_smf;

  int nsamples;
  int ncf;
  int iz;
  double satfac;
  double *mstar_threshold;
  double zlo, zhi;
  double mlo, mhi;
  double a[100];
  int *ndata;
  double *ngal;
  double *rdata[100];
  double *xdata[100];
  double *edata[100];
  double **covar[100];
  double **evect[100];
  double *eval[100];
  int npca;

  double mstar;
  double mstar_upper;
  double mstar_lower;
  int reset_inversion;
  int reset_fred;
  int mlow_flag;

  // this is for temp storage of HOD parameters
  double mmin, m1, mcut, mlow, slogm, alpha, maxcen;

  // this is for the red/blue samples
  double mstar_wplo[100][20];
  double mstar_wphi[100][20];

  int *ndatar;
  double *rdatar[100];
  double *xdatar[100];
  double *edatar[100];
  double **covarr[100];

  int *ndatab;
  double *rdatab[100];
  double *xdatab[100];
  double *edatab[100];
  double **covarb[100];


} wpl;


// this is for the lensing stuff (since many of the field 
// names will be the same as in the wpl structure
struct WPX {
  int nsamples;
  int ncf;
  int iz;
  double *mstar_threshold;
  double zlo, zhi;
  double mlo, mhi;
  double a[100];
  int *ndata;
  double *ngal;
  double *model[100];
  double *rdata[100];
  double *xdata[100];
  double *edata[100];
  double **covar[100];
  double mstar;
  double mstar_upper;
  double mstar_lower;
  int reset_inversion;
  int reset_fred;
  int mlow_flag;
  int calculate_two_halo;

  // this is for temp storage of HOD parameters
  double mmin, m1, mcut, mlow, slogm, alpha, maxcen;

  // this is for the red/blue samples
  double mstar_wplo[100][20];
  double mstar_wphi[100][20];

  int *ndatar;
  double *rdatar[100];
  double *xdatar[100];
  double *edatar[100];
  double **covarr[100];
  double *two_halo_red[100];

  int *ndatab;
  double *rdatab[100];
  double *xdatab[100];
  double *edatab[100];
  double **covarb[100];
  double *two_halo_blue[100];

} wpx;



/* Globals for LENSING analysis
 */
extern int SHMR_FLAG;             // flag for using SHMR rather than HOD
extern int LENSING_OUTPUT_FLAG;   // write_output_files
extern int HINVERSE_UNITS;        //
extern int COLOR;                 //
extern int SHMR_PARAMS;           //
extern int DONT_FIT_SMF;          // don't use the smf in the fit
extern int DONT_FIT_LENSING;      // don't use the lensing in the fit
extern int DONT_FIT_CLUSTERING;   // don't use the clustering in the fit
extern int BLUE_FLAG;             // are you modeling the blue subsample?
extern int JPL_FLAG;              // using the JPL cluster
extern int FISHER;                // using the Fisher matrix for covar (or to calculate it)
extern int RESET_FLAG_DS;         // reset flag for Delta Sigma
extern int RESET_FLAG_XGM1;       // reset flag for 1-halo xi_gm
extern int SATELLITE_PARAMETERIZATION; // flag for using the 5-parameter M1 function rather than the 2-parameter
extern int SMF_PCA;               // flag for using PCA analysis to get chi^2 for SMF model comparison
extern int VARIABLE_ALPHA;        // in SHMR analysis, if ALPHA can vary with stellar mass threshold
extern int VARIABLE_EXCLUSION;    // in SHMR analysis, have the exclusion radius be a variable 
extern int EXTERNAL_CONSTRAINTS;  // in mcmc_shmr, use extra data for constraining the functions 
extern int XIMATTER;              // flag for telling the Kaiser infall model to use matter clustering instead of galaxy
extern int NOKAISER;              // flag for not using RSD in wp calc
extern int CENTRALS_ONLY; 
