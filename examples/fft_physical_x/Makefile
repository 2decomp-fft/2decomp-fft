#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT) $(LFLAGS)

NP ?= 1
MPIRUN ?= mpirun

all: fft_c2c_x fft_r2c_x fft_grid_x

fft_c2c_x: fft_c2c_x.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

fft_r2c_x: fft_r2c_x.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

fft_grid_x: fft_grid_x.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_c2c_x
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_r2c_x
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_grid_x
else
check:
	$(MPIRUN) -n $(NP) ./fft_c2c_x
	$(MPIRUN) -n $(NP) ./fft_r2c_x
	$(MPIRUN) -n $(NP) ./fft_grid_x
endif

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all $(MPIRUN) -n 1 ./fft_c2c_x 1 1

clean:
	rm -f *.o fft_c2c_x fft_r2c_x fft_grid_x *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
