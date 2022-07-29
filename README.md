# 2decomp-fft
Customised 2DECOMP&amp;FFT library for Xcompact3d

Different compilers can be set by specifying `CMP`, e.g. `make CMP=intel`
to build with Intel compilers, see `Makefile` for options.

By default an optimised library will be built, debugging versions of the
library can be built with `make BUILD=debug`, a development version which 
additionally sets compile time flags to catch coding errors can be built 
with `make BUILD=dev` (GNU compilers only currently).
