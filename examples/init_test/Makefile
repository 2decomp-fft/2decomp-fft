FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LFLAGS)

OBJ = init_test.o

NP ?= 1
MPIRUN ?= mpirun

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
	rm -f *.o init_test *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
