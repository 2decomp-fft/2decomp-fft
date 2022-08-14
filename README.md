# 2decomp-fft

Different compilers can be set by specifying `CMP`, e.g. `make CMP=intel`
to build with Intel compilers, see `Makefile` for options.

By default an optimised library will be built, debugging versions of the
library can be built with `make BUILD=debug`, a development version which 
additionally sets compile time flags to catch coding errors can be built 
with `make BUILD=dev` (GNU compilers only currently).

## Testing and examples

Various example code to exercise 2decomp functionality can be found under ``examples/``
and can be built from the top-level directory by executing
```
make examples
```
which will (re)build 2decomp&fft as necessary.

**TODO** Convert examples to tests and automate running them
