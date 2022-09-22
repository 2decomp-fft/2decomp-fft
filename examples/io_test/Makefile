#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)

all: io_test io_read io_var_test io_plane_test io_bench

io_test: io_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_read: io_read.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_var_test: io_var_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_plane_test: io_plane_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_bench: io_bench.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

check:
	mpirun -n 12 ./io_test
	mpirun -n 6 ./io_read
	mpirun -n 6 ./io_var_test 3 2
	mpirun -n 12 ./io_plane_test
	mpirun -n 16 ./io_bench

clean:
	rm -f *.o io_test io_read io_var_test io_plane_test io_bench

realclean: clean
	rm -f *.dat io_var_data.* *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
