#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = timing.o

timing: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

check:
	mpirun -n 4 ./timing

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./timing 1 1

check:
	mpirun -n 4 ./timing

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./timing 1 1

clean:
	rm -f *.o timing *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
