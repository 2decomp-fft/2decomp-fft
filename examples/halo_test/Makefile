#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

OBJ = halo_test.o

halo_test: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 12 ./halo_test

clean:
	rm -f *.o halo_test

%.o : %.f90
	$(FC) $(INCLUDE) $(OPTIONS) $(FFLAGS) -c $<
