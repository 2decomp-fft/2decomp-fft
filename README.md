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
