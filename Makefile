#=======================================================================
# Makefile for 2DECOMP&FFT
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debugging
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

DEFS = -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
CMP = gcc# intel,gcc,nagfor,cray,nvhpc
FFT ?= generic# fftw3,fftw3_f03,generic,mkl,cufft
PARAMOD = mpi # multicore,gpu
PROFILER ?= none# none, caliper

BUILD ?= # debug can be used with gcc
FCFLAGS ?= # user can set default compiler flags
LDFLAGS ?= # user can set default linker flags
FFLAGS = $(FCFLAGS)
LFLAGS = $(LDFLAGS)
MODFLAG = -J 

LIBDECOMP = decomp2d

AR = ar
LIBOPT = rcs

#######CMP settings###########
CMPINC = Makefile.compilers
include $(CMPINC)

### List of files for the main code
SRCDECOMP = decomp_2d.f90 d2d_log.f90 io.f90

#######FFT settings##########
ifeq ($(FFT),fftw3)
  FFTW3_PATH ?= /usr
  FFTW3_PATH_INCLUDE ?= $(FFTW3_PATH)/include
  FFTW3_PATH_LIB ?= $(FFTW3_PATH)/lib
  INC=-I$(FFTW3_PATH_INCLUDE)
  LIBFFT=-L$(FFTW3_PATH_LIB) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH ?= /usr
  FFTW3_PATH_INCLUDE ?= $(FFTW3_PATH)/include
  FFTW3_PATH_LIB ?= $(FFTW3_PATH)/lib
  INC=-I$(FFTW3_PATH_INCLUDE)
  LIBFFT=-L$(FFTW3_PATH_LIB) -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  SRCDECOMP += ./glassman.f90
  INC=
  LIBFFT=
else ifeq ($(FFT),mkl)
  $(shell ls $(MKLROOT)/include && cp $(MKLROOT)/include/mkl_dfti.f90 src/)
  SRCDECOMP += mkl_dfti.f90
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC=-I$(MKLROOT)/include
else ifeq ($(FFT),cufft)
  CUFFT_PATH ?= $(NVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers
  INC=-I$(CUFFT_PATH)/include
  #LIBFFT=-L$(CUFFT_PATH)/lib64 -Mcudalib=cufft 
endif

### IO Options ###
LIBIO :=
OPTIO :=
INCIO :=
ADIOS2DIR :=
ifeq ($(IO),adios2)
  ifeq ($(ADIOS2DIR),)
    $(error Set ADIOS2DIR=/path/to/adios2/install/)
  endif
  OPTIO := -DADIOS2 $(OPT)
  INCIO := $(INC) $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-flags) #$(patsubst $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs),,$(shell $(ADIOS2DIR)/bin/adios2-config -f))
  LIBIO := $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs)
endif

### Add the profiler if needed
ifneq ($(PROFILER),none)
  DEFS += -DPROFILER
endif
ifeq ($(PROFILER),caliper)
  CALIPER_PATH ?= xxxxxxxxx/caliper/caliper_2.8.0
  SRCDECOMP := $(SRCDECOMP) profiler_caliper.f90
  INC := $(INC) -I$(CALIPER_PATH)/include/caliper/fortran
  LFLAGS := $(LFLAGS) -L$(CALIPER_PATH)/lib -lcaliper
endif

SRCDECOMP := $(SRCDECOMP) fft_$(FFT).f90
SRCDECOMP_ = $(patsubst %.f90,$(SRCDIR)/%.f90,$(SRCDECOMP))
OBJDECOMP = $(SRCDECOMP_:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

#######OPTIONS settings###########
OPT =
LINKOPT = $(FFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJDIR = obj
SRCDIR = src
DECOMPINC = mod
FFLAGS += $(MODFLAG)$(DECOMPINC) -I$(DECOMPINC)

OPT += $(OPTIO)
INC += $(INCIO)

include Makefile.settings

all: $(DECOMPINC) $(OBJDIR) $(LIBDECOMP)

$(DECOMPINC):
	mkdir $(DECOMPINC)

$(LIBDECOMP) : Makefile.settings lib$(LIBDECOMP).a

lib$(LIBDECOMP).a: $(OBJDECOMP)
	$(AR) $(LIBOPT) $@ $^

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDECOMP) : $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@

examples: $(LIBDECOMP)
	$(MAKE) -C examples

.PHONY: check

check: examples
	$(MAKE) -C examples $@

.PHONY: clean

clean: clean-examples
	rm -f $(OBJDIR)/*.o $(DECOMPINC)/*.mod $(DECOMPINC)/*.smod lib$(LIBDECOMP).a
	rm -f ./*.o ./*.mod ./*.smod # Ensure old files are removed
	rm -f Makefile.settings

clean-examples:
	$(MAKE) -C examples clean

.PHONY: Makefile.settings

Makefile.settings:
	echo "FC = $(FC)" > $@
	echo "FFLAGS = $(FFLAGS)" >> $@
	echo "OPT = $(OPT)" >> $@
	echo "DEFS = $(DEFS)" >> $@
	echo "INC = $(INC)" >> $@
	echo "LIBOPT = $(LIBOPT)" >> $@
	echo "LIBFFT = ${LIBFFT}" >> $@
	echo "LFLAGS = $(LFLAGS)" >> $@

export
