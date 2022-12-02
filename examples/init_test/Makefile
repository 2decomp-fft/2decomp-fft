FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LFLAGS)

OBJ = utilities.o init_test.o

NP ?= 1
MPIRUN ?= mpirun
MODDIR ?= $(DECOMPINC)

all: $(MODDIR) init_test

$(MODDIR):
	mkdir -p $@

init_test: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./init_test
else
check:
	$(MPIRUN) -n $(NP) ./init_test
endif

clean:
	rm -f *.o init_test *.log *.mod *.smod $(MODDIR)/*.mod $(MODDIR)/*.smod

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
