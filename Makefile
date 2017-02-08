
hd = $(HOME)/lib
LIB = -lm  #-lgsl -lgslcblas
#LIB = -lm -L$(HOME)/cosmo/lib -llatin #-lgsl -lgslcblas

#CC = /Library/Developer/CommandLineTools/usr/bin/gcc
CC = gcc
CFLAGS = -O2 
EXEC = HOD.x


# for ALEXIE
#OBJ_LENSING = chi2_lensing.o
OBJ_LENSING = chi2_lensing2.o one_halo_gm.o

OBJ_HOD = header.o main.o utility.o sigmac.o transfnc.o transfunc_file.o \
	nonlinear_power_spectrum.o least_squares.o test.o hod_functions.o \
	 xi_matter.o one_halo_rspace.o \
	input_params.o dFdx.o mstar.o halo_concentration.o growthfactor.o \
	halo_mass_conversion.o two_halo_rspace.o nfw_transform.o \
	 tasks.o wp_minimization.o \
	kaiser_distortions.o \
	tf_eisenstein_hu.o \
	populate_simulation_splash.o \
	dark_matter_statistics.o \
	meshlink2.o nbrsfind2.o i3tensor_2.o \
	wtheta.o sirko_integrate.o\
	populate_simulation_hod.o mcmc_shmr.o covar_pca.o shmr_minimization.o \
	external_constraints.o chain_postprocessing.o fit_for_bias.o \
	populate_simulation_hod_density2.o smf_primus.o

OBJ_CPP = latin_random.o

OBJ_EXTRA = sham.o convlv.o realft.o twofft.o four1.o
OBJ_NSF =  pvz_temp.o xi_multipoles.o one_halo_zspace.o two_halo_zspace.o galaxy_prob_vz.c \
	spherical_infall_velocity.o xi2d_interp.o

OBJ_SHMR = shmr_functions.o shmr_clustering.o 

OBJ_ERR = mcmc_with_errors.o halo_mass_function_error.o halo_bias_error.o
OBJ_STD = mcmc.o halo_mass_function.o halo_bias.o

OBJ_NR = nrutil.o qromo.o midpnt.o alexmidpnt.o alex1midpnt.o alex2midpnt.o \
	alex3midpnt.o midinf.o polint.o splint.o spline.o \
	zbrent.o qtrap.o trapzd.o alextrapzd.o cisi.o complex.o amoeba.o amotry.o \
	gaussj.o powell.o linmin.o f1dim.o mnbrak.o brent.o gasdev.o \
	ran1.o jacobi.o splin2.o splie2.o ran2.o sort2.o sort.o \
	svdcmp.o pythag.o sort2dbl.o

OBJ_LATIN = latin_hypercube.o

OBJS01 = $(OBJ_HOD) $(OBJ_NR) $(OBJ_STD) $(OBJ_LENSING) $(OBJ_SHMR) $(OBJ_EXTRA) $(OBJ_NSF) $(OBJ_LATIN)

$(EXEC): $(OBJS01)
#	gcc -c -o $(OBJS01) 
#	g++ -c -o latin_random.cpp latin_random.o
	$(CC) -o $@ $(OBJS01) $(LIB)
	cp -f $@ $(HOME)/exec/$@
$(OBJS01):	header.h nrutil.h complex.h


clean:
	rm -f *.o
