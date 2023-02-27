#include ../../src/Makefile.inc


#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LFLAGS)

NP ?= 1
MPIRUN ?= mpirun

all: test2d timing2d_real

test2d: test2d.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

timing2d_real: timing2d_real.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./test2d
	$(MPIRUN) -n $(NP) ./bind.sh ./timing2d_real
else
check:
	$(MPIRUN) -n $(NP) ./test2d
	$(MPIRUN) -n $(NP) ./timing2d_real
endif


clean:
	rm -f *.o test2d timing2d_real u*.dat *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
