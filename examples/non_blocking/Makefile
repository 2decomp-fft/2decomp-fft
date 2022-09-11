include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft $(LIBFFT)

all: blocking non_blocking

blocking: blocking.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

non_blocking: non_blocking.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o blocking non_blocking

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
