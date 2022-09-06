#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = timing.o

timing: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 4 ./timing

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./timing 1 1

clean:
	rm -f *.o timing

%.o : %.f90
	$(FC) $(INCLUDE) $(OPTIONS) $(FFLAGS) -c $<
