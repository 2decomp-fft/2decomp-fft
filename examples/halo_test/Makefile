#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

OBJ = halo_test.o

NP ?= 1
MPIRUN ?= mpirun

halo_test: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./halo_test
else
check:
	$(MPIRUN) -n $(NP) ./halo_test
endif

clean:
	rm -f *.o halo_test *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
