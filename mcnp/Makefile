# Fortran compiler (gnu or intel)
SYSTEM ?= gnu
DEBUG ?= false

LIBS = -llapack -lblas -lc -lm
LIBPATH =
INCLUDE =

# GNU Compiler Collection
ifeq ($(SYSTEM),gnu)
    FC = gfortran
    FFLAGS = -Wall -Wsurprising -Wconversion -Wline-truncation -std=f2008 -Wno-maybe-uninitialized
    FFLAGS_DEBUG = -g -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow,underflow -fimplicit-none -Wcharacter-truncation -Waliasing -Wimplicit-interface  -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -fcheck=all -fopenmp
    FFLAGS_OPT = -O3 -march=native -fwhole-file -fopenmp
    F77FLAGS = -g -gdwarf-2 -fbacktrace
    F77FLAGS_DEBUG = -finit-real=snan -ffpe-trap=invalid,overflow,zero
    F77FLAGS_OPT = -O3 -funroll-loops -ftree-vectorize
endif

# Intel Compiler Suite
ifeq ($(SYSTEM),intel)
    FC = ifx
    FFLAGS = -stand f08 -warn all -qopenmp -heap-arrays 1024
    FFLAGS_DEBUG = -O0 -g -check all -debug extended -warn all -fpe0 -traceback -ftrapuv -qmkl=sequential
    FFLAGS_OPT = -Ofast -xHost -fno-alias -parallel -qmkl=sequential
    F77FLAGS = -g -traceback -fpe0
    F77FLAGS_DEBUG = -O0 -ftrapuv
    F77FLAGS_OPT = -O2 -xHost
    LIBPATH += -L$(MKLROOT)/lib
    INCLUDE = -I$(MKLROOT)/include
    LIBS = -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

# Debug flags
ifeq ($(DEBUG),true)
    FFLAGS += $(FFLAGS_DEBUG)
    F77FLAGS += $(F77FLAGS_DEBUG)
else
    FFLAGS += $(FFLAGS_OPT)
    F77FLAGS += $(F77FLAGS_OPT)
endif

LBFGSB_OBJS = \
	lbfgsb/lbfgsb.o \
	lbfgsb/linpack.o \
	lbfgsb/timer.o

MODEL_OBJS = \
	osl.o \
	encoding.o \
	stdlib_experimental_sparse.o \
	sparse.o

MCNP_OBJS = \
	osl.o \
	encoding.o \
	stdlib_experimental_sparse.o \
	sparse.o \
	dataset.o \
	model.o \
	lbfgsb.o \
	$(LBFGSB_OBJS)

TARGETS = mcnp

all: $(TARGETS)

encoding.o: encoding.f90
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

osl.o: osl.f90
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

stdlib_experimental_sparse.o: stdlib_experimental_sparse.f90 osl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

sparse.o: sparse.f90 osl.o stdlib_experimental_sparse.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

dataset.o: dataset.f90 osl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

model.o: model.f90 $(MODEL_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

lbfgsb.o: lbfgsb.f90 osl.o
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/lbfgsb.o: lbfgsb/lbfgsb.f
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/linpack.o: lbfgsb/linpack.f
	$(FC) $(F77FLAGS) -o $@ -c $<

lbfgsb/timer.o: lbfgsb/timer.f
	$(FC) $(F77FLAGS) -o $@ -c $<

mcnp.o: mcnp.f90 osl.o model.o lbfgsb.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

mcnp: $(MCNP_OBJS) mcnp.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(LIBPATH) $(MCNP_OBJS) mcnp.o $(LIBS)

clean:
	rm -f *.mod *.dat *.o $(TARGETS) $(TARGETS:=.o)
	rm -f $(MODEL_OBJS) $(DATASET_OBJS) $(MCNP_OBJS)
	rm -f $(LBFGSB_OBJS)
	-rm -rf $(TARGETS:=.dSYM)
	-rm -rf *__genmod.f90
