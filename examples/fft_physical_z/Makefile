#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT) $(LFLAGS)

NP ?= 1
MPIRUN ?= mpirun

all: fft_c2c_z fft_r2c_z 

fft_c2c_z: fft_c2c_z.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

fft_r2c_z: fft_r2c_z.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_c2c_z
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_r2c_z
else
check:
	$(MPIRUN) -n $(NP) ./fft_c2c_z
	$(MPIRUN) -n $(NP) ./fft_r2c_z
endif

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all $(MPIRUN) -n 1 ./fft_c2c_z 1 1

clean:
	rm -f *.o fft_c2c_z fft_r2c_z *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
