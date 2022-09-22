#include ../../src/Makefile.inc


#INCLUDE = -I../../include
blank :=
space := $(blank) $(blank)
MODFLAG_NSP := $(subst $(space),,$(MODFLAG))
FFLAGS := $(subst ?,$(space),$(subst $(MODFLAG_NSP)?,$(MODFLAG_NSP)?../../,$(subst $(space),?,$(FFLAGS))))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

OBJ = test2d.o

NP ?= 1
MPIRUN ?= mpirun

test2d: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	$(MPIRUN) -n $(NP) ./test2d

clean:
	rm -f *.o test2d u*.dat *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
