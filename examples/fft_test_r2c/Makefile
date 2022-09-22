#include ../../src/Makefile.inc

#INCLUDE = -I../../include
blank :=
space := $(blank) $(blank)
MODFLAG_NSP := $(subst $(space),,$(MODFLAG))
FFLAGS := $(subst ?,$(space),$(subst $(MODFLAG_NSP)?,$(MODFLAG_NSP)?../../,$(subst $(space),?,$(FFLAGS))))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = fft_test_r2c.o

NP ?= 1
MPIRUN ?= mpirun

fft_test_r2c: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	$(MPIRUN) -n $(NP) ./fft_test_r2c

clean:
	rm -f *.o fft_test_r2c *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
