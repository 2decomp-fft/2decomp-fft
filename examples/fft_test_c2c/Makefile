#include ../../src/Makefile.inc

#INCLUDE = -I../../include
blank :=
space := $(blank) $(blank)
FFLAGS := $(subst ?,$(space),$(subst -J?,-J?../../,$(subst $(space),?,$(FFLAGS))))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = fft_test_c2c.o

MPIRUN ?= mpirun
NP ?= 1

fft_test_c2c: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	$(MPIRUN) -n $(NP) ./fft_test_c2c

clean:
	rm -f *.o fft_test_c2c *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $^ -o $@
