#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = fft_test_c2c.o

fft_test_c2c: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o fft_test_c2c

%.o : %.f90
	$(FC) $(INC) $(OPTIONS) $(FFLAGS) -c $<
