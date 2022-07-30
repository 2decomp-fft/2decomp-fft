#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J%,-J../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = timing.o

timing: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o timing

%.o : %.f90
	$(FC) $(INCLUDE) $(OPTIONS) $(FFLAGS) -c $<
