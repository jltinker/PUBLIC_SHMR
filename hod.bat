%  Cosmological Parameters
%----------------------------------------------------------------------------
%
GAMMA		0.55
OMEGA_M		0.27
SIGMA_8		0.82
RHO_CRIT	2.775E11		% h^2 M_sol/Mpc^3
SPECTRAL_INDX	0.95
HUBBLE		0.7
OMEGA_B		0.04
DELTA_CRIT	1.686
ITRANS		5
TF_file		transfunc.WMAP3
REDSHIFT	0

LENSING_OUTPUT_FLAG	0
STEPFAC			0.6
FISHER			5
JPL_FLAG		0
DONT_FIT_LENSING	0
DONT_FIT_CLUSTERING	0
DONT_FIT_SMF		0

% N-body simulation specifics
%----------------------------------------------------------------------------

DELTA_HALO	360
BOX_SIZE	250		% Mpc/h
RESOLUTION	0.375		% BOX_SIZE/Npart^(1/3)


% HOD Parameters
%----------------------------------------------------------------------------
% NB! -> If M_min is > 0, then GALAXY_DENSITY is calculated in the code.
%	 If M_min is <=0, then M_min is calculated for the given GALAXY_DENSITY

GALAXY_DENSITY  0.00308     

M_min           0
M1              2.91e13
M_max           1.00E+15
M_cut           1.087e7
alpha           1.20
sigma_logM 	0.15
pdfs		11
pdfc		2

VBIAS		1
VBIAS_C		0
CVIR_FAC	1
EXCLUSION	4

% Free parameters for HOD fitting
% -----------------------------------------------------------------------------

free[1]		0	% M_min
free[2]		1	% M1
free[3]		1	% alpha
free[4]		1	% M_cut
free[5]		0	% sigma_logM
free[6]		0	% CVIR_FAC
free[7]		0	% HOD.MaxCen


% Parameters for the Stellar-to-Halo Mass Relation
% ------------------------------------------------------------------------------

SHMR_FLAG	1	% Use the SHMR rather than the HOD 
COLOR		0	% Use color-dependent SHMR (red and blue)

SHMR.MhaloNorm	12.73
SHMR.MstarNorm	11.04
SHMR.beta	0.466
SHMR.delta	0.61
SHMR.gamma	1.95
SHMR.scatter	0.25
SHMR.McutNorm	1.65
SHMR.M1Norm	9.04
SHMR.McutSlope	0.59
SHMR.M1Slope	0.74
SHMR.alpha_sat	1

SHMR2.MhaloNorm	12.62
SHMR2.MstarNorm	10.92
SHMR2.beta	0.457
SHMR2.delta	0.766
SHMR2.gamma	1.53
SHMR2.scatter	0.206
SHMR2.McutNorm	1.47
SHMR2.M1Norm	20.62
SHMR2.McutSlope	-0.13
SHMR2.M1Slope	0.859
SHMR2.alpha_sat	1

red_central_fraction[1] -7.088729e+00 
red_central_fraction[2] -1.120370e+00 
red_central_fraction[3] 4.136595e-01 
red_central_fraction[4] 6.680107e-01 
red_central_fraction[5] 7.951682e-01 

% If you are fitting the SHMR to data, please specify the files that contain
% the names of the files that contain the data, as well as the file that contains
% 
% ------------------------------------------------------------------------------

StellarMassBinsClustering	clustering_bins.dat
StellarMassBinsLensing		clustering_bins.dat
ClusteringFilenames		xxxxx
StellarMassFunctionFilename	xxxxx
LensingFilenames		xxxxx


% These are the different tasks that the program can perform
% If the All flag is hit, then the other flags are ignored and all the tasks
% are performed.
%-------------------------------------------------------------------------------

clustering	0
lensing		1
wp_minimize	0
HOD		0
root_filename	test_z0
populate_sim	0
SHMR 		0

halostats	0
matter_xi	0
matter_pk	0
sigma_r		0
cvir		0

HaloFile	/Users/tinker/cosmo/HALO_CATALOGS/BOLSHOI/halo_bolshoi_z0.00_D360.dat

fname_wp	/home/tinker/cosmo/SDSS/Data/wp_btm20.5L.dat
fname_covar	/home/tinker/cosmo/SDSS/Data/wp_btm20.5L_covar.dat

% Files and flags for zspace data and minimization
%-------------------------------------------------------------------------------

POWELL		0
MCMC		0
COVAR		0

OUTPUT		1

