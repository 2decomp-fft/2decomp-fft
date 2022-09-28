#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT) $(LFLAGS)

OBJ = fft_physical_x.o

fft_physical_x: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 4 ./fft_physical_x

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./fft_physical_x 1 1

check:
	$(MPIRUN) -n $(NP) ./fft_test_r2c

clean:
	rm -f *.o fft_physical_x *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
