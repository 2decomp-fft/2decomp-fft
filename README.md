# 2decomp-fft

## Building

Different compilers can be set by specifying `CMP`, e.g. `make CMP=intel`
to build with Intel compilers, see `Makefile` for options.

By default an optimised library will be built, debugging versions of the
library can be built with `make BUILD=debug`, a development version which 
additionally sets compile time flags to catch coding errors can be built 
with `make BUILD=dev` (GNU compilers only currently). The behavior of debug
and development versions of the library can be changed before the initialization
using the variable ``decomp_debug`` or the environment variable ``DECOMP_2D_DEBUG``.
The value provided with the environment variable must be a positive integer below 9999.

On each build of the library (`make`, `make all`) a temporary file `Makefile.settings` with
all current options (`FFLAGS`, `DEFS`, etc.) will be created, and included
on subsequent invocations, the user therefore does not need to keep
specifying options between builds.

To perform a clean build run `make clean` first, this will delete all
output files, including `Makefile.settings`.

## Testing and examples

Various example code to exercise 2decomp functionality can be found under ``examples/``
and can be built from the top-level directory by executing
```
make check
```
which will (re)build 2decomp&fft as necessary.

**TODO** Convert examples to tests and automate running them

## Profiling

Profiling can be activated in the Makefile. Set the variable `PROFILER` to one of the supported profilers (only `caliper` currently). If using `caliper`, provide the installation path in the variable `CALIPER_PATH`. When the profiling is active, one can tune it before calling `decomp_2d_init` using the subroutine `decomp_profiler_prep`. The input argument for this subroutine is a logical array of size 4. Each input allow activation / deactivation of the profiling as follows :

1. Profile transpose operations (default : true)
2. Profile IO operations (default : true))
3. Profile FFT operations (default : true)
4. Profile decomp_2d init / fin subroutines (default : true)

## Optional dependencies

### FFTW

TODO

### Caliper

The library [caliper](https://github.com/LLNL/Caliper) can be used to profile the execution of the code. The version 2.8.0 was tested and is supported. Please note that one must build caliper and decomp2d against the same C/C++/Fortran compilers and MPI libray. For build instructions, please check [here](https://github.com/LLNL/Caliper#building-and-installing) and [here](https://software.llnl.gov/Caliper/CaliperBasics.html#build-and-install). Below is a suggestion for the compilation of the library using the GNU compilers:

```
git clone https://github.com/LLNL/Caliper.git caliper_github
cd caliper_github
git checkout v2.8.0
mkdir build && cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../caliper_build_2.8.0 -DWITH_FORTRAN=yes -DWITH_MPI=yes -DBUILD_TESTING=yes ../
make -j
make test
make install
```
