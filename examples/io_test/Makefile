#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(patsubst -J\ %,-J\ ../../%,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP)
ifneq (,$(findstring DT3PIO,$(OPTIONS)))
  LIBS+= -L$(T3PIO_PATH)/lib -lt3pio
endif

all: io_test io_read io_var_test io_plane_test io_bench

io_test: io_test.o
	$(FC) -o $@ $@.o $(LIBS)

io_read: io_read.o
	$(FC) -o $@ $@.o $(LIBS)

io_var_test: io_var_test.o
	$(FC) -o $@ $@.o $(LIBS)

io_plane_test: io_plane_test.o
	$(FC) -o $@ $@.o $(LIBS)

io_bench: io_bench.o
	$(FC) -o $@ $@.o $(LIBS)

clean:
	rm -f *.o io_test io_read io_var_test io_plane_test io_bench

realclean: clean
	rm -f *.dat io_var_data.*

%.o : %.f90
	$(FC) $(INCLUDE) $(OPTIONS) $(ARG) $(FFLAGS) -c $<
