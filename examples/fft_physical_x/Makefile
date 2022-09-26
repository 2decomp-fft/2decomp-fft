#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LIBFFT)

OBJ = fft_physical_x.o

fft_physical_x: $(OBJ)
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	mpirun -n 4 ./bind.sh ./fft_physical_x
else
check:
	mpirun -n 4 ./fft_physical_x
endif

mem_leak:
	valgrind --leak-check=full --show-leak-kinds=all mpirun -n 1 ./fft_physical_x 1 1

clean:
	rm -f *.o fft_physical_x *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
