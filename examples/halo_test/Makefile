#include ../../src/Makefile.inc

#INCLUDE = -I../../include
blank :=
space := $(blank) $(blank)
FFLAGS := $(subst ?,$(space),$(subst -J?,-J?../../,$(subst $(space),?,$(FFLAGS))))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

OBJ = halo_test.o

NP ?= 1
MPIRUN ?= mpirun

halo_test: $(DECOMPINC) $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

$(DECOMPINC):
	mkdir $@

check:
	$(MPIRUN) -n $(NP) ./halo_test

clean:
	rm -f *.o halo_test *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
