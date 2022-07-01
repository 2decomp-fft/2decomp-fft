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
FFT = generic# fftw3,fftw3_f03,generic,mkl
PARAMOD = mpi # multicore,gpu

BUILD ?= # debug can be used with gcc
FCFLAGS ?= # user can set default compiler flags
LDFLAGS ?= # user can set default linker flags
FFLAGS = $(FCFLAGS)
LFLAGS = $(LDFLAGS)

LIBDECOMP = libdecomp2d.a

AR = ar
LIBOPT = rcs

#######CMP settings###########
CMPINC = Makefile.compilers
include $(CMPINC)

### List of files for the main code
SRCDECOMP = ./decomp_2d.f90 ./glassman.f90 ./fft_$(FFT).f90 
OBJDECOMP = $(SRCDECOMP:%.f90=$(OBJDIR)/%.o)

#######FFT settings##########
ifeq ($(FFT),fftw3)
  FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH=/usr                                #ubuntu # apt install libfftw3-dev
  #FFTW3_PATH=/usr/lib64                         #fedora # dnf install fftw fftw-devel
  #FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1     #macOS  # brew install fftw
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  INC=
  LIBFFT=
else ifeq ($(FFT),mkl)
  SRCDECOMP := $(DECOMPDIR)/mkl_dfti.f90 $(SRCDECOMP)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC=-I$(MKLROOT)/include
else ifeq ($(FFT),cufft)
  #CUFFT_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.1/math_libs                                
  INC=-I${NVHPC}/Linux_x86_64/${EBVERSIONNVHPC}/compilers/include
  #LIBFFT=-L$(CUFFT_PATH)/lib64 -Mcudalib=cufft 
endif

#######OPTIONS settings###########
OPT =
LINKOPT = $(FFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJDIR = obj
DECOMPINC = mod
FFLAGS += -J$(DECOMPINC) -I$(DECOMPINC)

all: $(DECOMPINC) $(OBJDIR) $(LIBDECOMP)

$(DECOMPINC):
	mkdir $(DECOMPINC)

$(LIBDECOMP) : $(OBJDECOMP)
	$(AR) $(LIBOPT) $@ $^

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDECOMP) : $(OBJDIR)/%.o : ./%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJDECOMP) $(DECOMPINC)/*.mod $(LIBDECOMP)
