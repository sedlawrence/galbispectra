# Makefile for Song.
# Last edited by Guido Walter Pettinari on 27/10/2013.


# =========================================================================
# =                         Compilers and flags                           =
# =========================================================================

# Your C compiler and library tool
CC			 = /usr/local/bin/gcc -g
AR       = ar rv

# Your Fortran compiler
FC       = /usr/local/bin/gfortran

# Optimization flags.
# The flags -fast, -O3 and -ffast-math seem to have the same effect, that is a factor 2 speed-up.
# The fastest combination seems to be -O3 -ffast-math, whic delivers a speed-up of ~3.
# Note that (i) the flag -ffast-math will produce different results at the 1 permille level,
# and that (ii) debug won't work if the optimization is on.
OPTFLAG         = -O
FORTRAN_OPTFLAG = -O

# Openmp flag (comment for compiling without openmp)
CCFLAG     = -fopenmp -std=c99



# =========================================================================
# =                        Directories locations                          =
# =========================================================================


# Source files
vpath %.c source:tools:main:test
vpath %.f90 source:tools:main:test
vpath %.f source:tools:main:test
vpath %.o build
vpath .base build

# Build directory
MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

# Header files and libraries.
# Remember to include gfortran or ifort libraries in order to use Slatec.
INCLUDES 					  = -I../include
FORTRAN_LIB_DIR 		=  
LIBRARIES           = -lgfortran -fopenmp -lm







# =========================================================================
# =                          Compilation rules                            =
# =========================================================================

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o:  %.f90 .base
	cd $(WRKDIR);$(FC) $(FFLAGS) $(FORTRAN_OPTFLAG) -c ../$< -o $*.o

%.o:  %.f .base
	cd $(WRKDIR);$(FC) $(FFLAGS) $(FORTRAN_OPTFLAG) -c ../$< -o $*.o



# ==========================================================================
# =                             Object files                               =
# ==========================================================================

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o evolver_ndf15.o arrays.o parser.o quadrature.o
INPUT = input.o
INPUT2 = input2.o
PRECISION = precision.o
BACKGROUND = background.o
THERMO = thermodynamics.o
PERTURBATIONS = perturbations.o
PERTURBATIONS2 = perturbations2.o
BESSEL = bessel.o
BESSEL2 = bessel2.o
TRANSFER = transfer.o
TRANSFER2 = transfer2.o
PRIMORDIAL = primordial.o
SPECTRA = spectra.o
BISPECTRA = bispectra.o
BISPECTRA2 = bispectra2.o
FISHER = fisher.o
NONLINEAR = trg.o nonlinear.o
LENSING = lensing.o
OUTPUT = output.o
SONG_TOOLS = $(TOOLS) song_tools.o slatec_3j_f90.o mesh_interpolation.o




# ==========================================================================
# =                            Executables                                 =
# ==========================================================================

# Rule to be executed when make is called with no parameters
default: class song print_params print_sources1 print_sources2 print_transfers1 print_transfers2 print_bispectra

# CLASS code
CLASS = class.o
TEST_PLEGENDRE = test_plegendre.o
TEST_LOOPS = test_loops.o
TEST_TRANSFER = test_transfer.o
TEST_PERTURBATIONS = test_perturbations.o
TEST_THERMODYNAMICS = test_thermodynamics.o
TEST_BACKGROUND = test_background.o

libclass.a: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT)
	$(AR)  $@ $(addprefix build/,$(notdir $^))
 
class: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) $(FISHER) $(NONLINEAR) $(LENSING) $(OUTPUT) $(CLASS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_plegendre: $(TOOLS) $(SONG_TOOLS) $(TEST_PLEGENDRE)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(TEST_LOOPS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(TEST_TRANSFER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(TEST_PERTURBATIONS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(TEST_THERMODYNAMICS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(INPUT) $(BACKGROUND) $(TEST_BACKGROUND)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

clean: .base
	rm -rf $(WRKDIR);


# SONG code
SONG = song.o
PRINT_K = print_k.o
PRINT_K_2ND_ORDER = print_k_2nd_order.o
PRINT_TAU_2ND_ORDER = print_tau_2nd_order.o
PRINT_PARAMS = print_params.o
PRINT_SOURCES1 = print_sources1.o
PRINT_SOURCES2 = print_sources2.o
PRINT_TRANSFERS1 = print_transfers1.o
PRINT_TRANSFERS2 = print_transfers2.o
PRINT_BISPECTRA = print_bispectra.o
PRINT_BACKGROUND = print_background.o
PRINT_THERMO = print_thermo.o
PRINT_INITIAL_CONDITIONS_2ND_ORDER = print_initial_conditions_2nd_order.o
PRINT_KERNEL = print_kernel.o
PRINT_MATRICES = print_matrices.o
TEST_COUPLINGS = test_couplings.o

test_bispectra: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) test_bispectra.o
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

song: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) $(BISPECTRA2)\
	$(FISHER) $(NONLINEAR) $(LENSING) $(OUTPUT) $(SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_params: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_PARAMS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources1: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2) $(TRANSFER) $(TRANSFER2) $(PRINT_SOURCES1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources2: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(PRINT_SOURCES2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers1: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRINT_TRANSFERS1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers2: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRINT_TRANSFERS2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_bispectra: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) $(BISPECTRA2) $(FISHER) $(PRINT_BISPECTRA)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_background: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PRINT_BACKGROUND)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_thermo: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PRINT_THERMO)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_K)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_K_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_tau_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_TAU_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_initial_conditions_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_INITIAL_CONDITIONS_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

print_kernel: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_KERNEL)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_matrices: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THREEJ) $(PRINT_MATRICES)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

test_couplings: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(TEST_COUPLINGS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm