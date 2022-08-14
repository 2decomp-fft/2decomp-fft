#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

OBJ = test2d.o

test2d: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 12 ./test2d

clean:
	rm -f *.o test2d u*.dat

%.o : %.f90
	$(FC) $(INC) $(OPTIONS) $(FFLAGS) -c $<
