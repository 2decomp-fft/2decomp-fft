#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT) $(LFLAGS)

OBJ = fft_physical_z.o

NP ?= 4
MPIRUN ?= mpirun

fft_physical_z: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./fft_physical_z
else
check:
	$(MPIRUN) -n $(NP) ./fft_physical_z
endif

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all $(MPIRUN) -n 1 ./fft_physical_z 1 1

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./timing 1 1

clean:
	rm -f *.o fft_physical_z *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
