#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT) $(LFLAGS)

OBJ = fft_physical_z.o

fft_physical_z: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 4 ./fft_physical_z

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./fft_physical_z 1 1

check:
	mpirun -n 4 ./timing

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./timing 1 1

clean:
	rm -f *.o fft_physical_z *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
