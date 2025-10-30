# Fortran compiler (gnu or intel)
SYSTEM ?= gnu
DEBUG ?= false

# GNU Compiler Collection
ifeq ($(SYSTEM),gnu)
    FC = gfortran
    FFLAGS = -g -gdwarf-2 -fbacktrace -Wall -Wsurprising -Wconversion -Wno-unused -std=f2008
    FFLAGS_DEBUG = -fbounds-check -finit-real=snan -ffpe-trap=invalid,overflow,zero -fopenmp
    FFLAGS_OPT = -O3 -funroll-loops -ftree-vectorize -fopenmp
    F77FLAGS = -g -gdwarf-2 -fbacktrace
    F77FLAGS_DEBUG = -finit-real=snan -ffpe-trap=invalid,overflow,zero
    F77FLAGS_OPT = -O3 -funroll-loops -ftree-vectorize
endif

# Intel Compiler Suite
ifeq ($(SYSTEM),intel)
    FC = ifort
    FFLAGS = -stand f08 -g -traceback -fpe0 -qopenmp
    FFLAGS_DEBUG = -O0 -check all -debug extended -warn all -ftrapuv
    FFLAGS_OPT = -O2 -xHost
    F77FLAGS = -g -traceback -fpe0
    F77FLAGS_DEBUG = -O0 -ftrapuv
    F77FLAGS_OPT = -O2 -xHost
endif

# Intel Math Kernel Library or standard LAPACK and BLAS
ifeq ($(MKL),true)
    FFLAGS += -mkl
    F77FLAGS += -mkl
else
    LIBS = -llapack -lblas -lm
endif

# Debug flags
ifeq ($(DEBUG),true)
    FFLAGS += $(FFLAGS_DEBUG)
    F77FLAGS += $(F77FLAGS_DEBUG)
else
    FFLAGS += $(FFLAGS_OPT)
    F77FLAGS += $(F77FLAGS_OPT)
endif

INCLUDE = -I.

MODEL_OBJS = \
	osl.o \
	dataset.o \
	rust_data.o \
	matrix.o

MC1P_OBJS = \
	osl.o \
	expokit/dgpadm.o \
	lbfgsb/timer.o \
	lbfgsb/linpack.o \
	lbfgsb/lbfgsb.o \
	lbfgsb.o \
	dataset.o \
	matrix.o \
	rust_data.o \
	rust_model.o

RUSTCT_OBJS = $(MC1P_OBJS)

TARGETS = mc1p rustct-heterogeneous rustct-homogeneous rustct-abbe

all: $(TARGETS)

osl.o: osl.f90
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

expokit/%.o: expokit/%.f
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/linpack.o: lbfgsb/linpack.f
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/timer.o: lbfgsb/timer.f
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/lbfgsb.o: lbfgsb/lbfgsb.f lbfgsb/timer.o lbfgsb/linpack.o
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb.o: lbfgsb.f90 lbfgsb/lbfgsb.o osl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

matrix.o: matrix.f90 osl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

rust_data.o: rust_data.f90 dataset.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

rust_model.o: rust_model.f90 $(MODEL_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

dataset.o: dataset.f90 osl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

mc1p.o: mc1p.f90 $(MC1P_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

mc1p: $(MC1P_OBJS) mc1p.o
	$(FC) $(FFLAGS) -o $@ $(LIBPATH) $(MC1P_OBJS) mc1p.o $(LIBS)

rustct-heterogeneous: $(RUSTCT_OBJS) rustct-heterogeneous.o
	$(FC) $(FFLAGS) -o $@ $(LIBPATH) $(RUSTCT_OBJS) rustct-heterogeneous.o $(LIBS)

rustct-heterogeneous.o: rustct.f90 $(RUSTCT_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -cpp -o $@ -c $<

rustct-homogeneous: $(RUSTCT_OBJS) rustct-homogeneous.o
	$(FC) $(FFLAGS) -o $@ $(LIBPATH) $(RUSTCT_OBJS) rustct-homogeneous.o $(LIBS)

rustct-homogeneous.o: rustct.f90 $(RUSTCT_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -cpp -DMODEL_HOMOGENEOUS -o $@ -c $<

rustct-abbe: $(RUSTCT_OBJS) rustct-abbe.o
	$(FC) $(FFLAGS) -o $@ $(LIBPATH) $(RUSTCT_OBJS) rustct-abbe.o $(LIBS)

rustct-abbe.o: rustct.f90 $(RUSTCT_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -cpp -DMODEL_ABBE -o $@ -c $<

clean:
	rm -f *.mod *.dat $(TARGETS) $(TARGETS:=.o)
	rm -f rustct-heterogeneous rustct-homogeneous rustct-abbe
	rm -f rustct-heterogeneous.o rustct-homogeneous.o rustct-abbe.o
	rm -f $(MODEL_OBJS) $(DATASET_OBJS) $(MC1P_OBJS) $(RUSTCT_OBJS)
	-rm -rf $(TARGETS:=.dSYM) rustct-heterogeneous.dSYM rustct-homogeneous.dSYM rustct-abbe.dSYM
