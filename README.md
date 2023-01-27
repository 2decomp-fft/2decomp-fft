# 2decomp-fft

## Building

The build system is driven by `cmake`, to generate a build system run
```
cmake -B build .
```
from the `2decomp&fft` root directory.
It may be necessary to override the Fortran compiler if the wrong one is picked up (due to, e.g.
`PATH` ordering), for example to specify the `MPICH` compiler on Ubuntu
```
FC=mpif90.mpich cmake -B build .
```
This should generate a directory `build/` containing the build system, this can be further
configured using the `ccmake` utility
```
ccmake build/
```
and editing as desired, variables that are likely of interest are: `CMAKE_BUILD_TYPE` and `FFT_Choice`;
additional variables can be shown by entering "advanced mode" by pressing `t`.
By default a `RELEASE` build will built, other options for `CMAKE_BUILD_TYPE` are `DEBUG` and `DEV` which
turn on debugging flags and additionally try to catch coding errors at compile time, respectively.
The behavior of debug and development versions of the library can be changed before the
initialization using the variable ``decomp_debug`` or the environment variable ``DECOMP_2D_DEBUG``.
The value provided with the environment variable must be a positive integer below 9999.

Once the build system has been configured, build `2decomp&fft` by running
```
make -j <nproc> -C build
```
appending ```VERBOSE=1``` will display additional information about the build, such as compiler flags.
After building the library can be tested by running
```
make -C build test
```
which uses the `ctest` utility.
Finally, install the library by running
```
make -C build install
```
which will install the library under `build/opt` by default, this location can be configured using
`ccmake` as described above and setting the variable `CMAKE_INSTALL_PREFIX`.

Occasionally a clean build is required, this can be performed by running
```
make -C build clean
```

## Testing and examples

Various example code to exercise 2decomp functionality can be found under ``examples/``
and can be built from the top-level directory by executing
```
make check
```
which will (re)build 2decomp&fft as necessary.

By default this will run the tests/examples on a single rank, launched with ``mpirun``.
The number of ranks can be controlled by setting the ``NP`` variable and the launcher is overridden by
setting ``MPIRUN``, for example to run tests on 4 ranks using an ``mpiexec`` launcher installed at
``${HOME}/bin/mpiexec`` the command would be
```
make NP=4 MPIRUN=${HOME}/bin/mpiexec check
```

## GPU compilation

The library can perform multi GPU offoloading using the NVHPC compiler suite for NVIDIA hardware. 
The implementation is based on CUDA-aware MPI and NVIDIA Collective Communication Library (NCCL).
The FFT is based on cuFFT. 
To compile the library for GPU it is possible to execute the following
```
make CMP=nvhpc FFT=cufft PARAMOD=gpu CUFFT_PATH=PATH_TO_NVHPC/Vers/Linux_x86_64/Vers/compilers/ 
``` 
The `Makefile` will look for the relative libraries (NVCC, cuFFT, etc) under the `${CUFFT_PATH}/include`
NCCL is not activated by default. If NCCL is installed/required use `NCCL=yes`. 
The current implementation relays also on opeanACC
and on automatic optimization of `do concurrent` loops.
By default the compute architecture for the GPU is 80 (i.e. Ampere), to change it use `CCXY=XY` 
 
## Profiling

Profiling can be activated in the Makefile. Set the variable `PROFILER` to one of the supported profilers (only `caliper` currently). If using `caliper`, provide the installation path in the variable `CALIPER_PATH`. When the profiling is active, one can tune it before calling `decomp_2d_init` using the subroutine `decomp_profiler_prep`. The input argument for this subroutine is a logical array of size 4. Each input allow activation / deactivation of the profiling as follows :

1. Profile transpose operations (default : true)
2. Profile IO operations (default : true)
3. Profile FFT operations (default : true)
4. Profile decomp_2d init / fin subroutines (default : true)

## FFT backends

The library provides a built-in FFT engine and supports various FFT backends : FFTW, Intel oneMKL, Nvidia cuFFT. The FFT engine selected during compilation is available through the variable `D2D_FFT_BACKEND` defined in the module `decomp_2d_fft`. The expected value is defined by the integer constants
```
integer, parameter, public :: D2D_FFT_BACKEND_GENERIC = 0   ! Built-in engine
integer, parameter, public :: D2D_FFT_BACKEND_FFTW3 = 1     ! FFTW
integer, parameter, public :: D2D_FFT_BACKEND_FFTW3_F03 = 2 ! FFTW (Fortran 2003)
integer, parameter, public :: D2D_FFT_BACKEND_MKL = 3       ! Intel oneMKL
integer, parameter, public :: D2D_FFT_BACKEND_CUFFT = 4     ! Nvidia cuFFT
```
exported by the module `decomp_2d_constants`.
The external code can use the named variables to check the FFT backend used in a given build.

### OVERWRITE flag

- The generic backend supports the OVERWRITE flag but it can not perform in-place transforms
- The FFTW3 and FFTW3_F03 backends support the OVERWRITE flag and can perform in-place complex 1D fft
- The oneMKL backend supports the OVERWRITE flag and can perform in-place complex 1D fft
- The cuFFT backend supports the OVERWRITE flag and can perform in-place complex 1D fft

## Miscellaneous

### Print the log to a file or to stdout

Before calling `decomp_2d_init`, the external code can modify the variable `decomp_log` to change the output for the log. The expected value is defined by the integer constants
```
integer, parameter, public :: D2D_LOG_QUIET = 0       ! No logging output
integer, parameter, public :: D2D_LOG_STDOUT = 1      ! Root rank logs output to stdout
integer, parameter, public :: D2D_LOG_TOFILE = 2      ! Root rank logs output to the file "decomp_2d_setup.log"
integer, parameter, public :: D2D_LOG_TOFILE_FULL = 3 ! All ranks log output to a dedicated file
```
exported by the `decomp_2d_constants` module.
Although their values are shown here, users should not rely on these and are recommended to prefer to use the named variables `D2D_LOG_QUIET`, etc. instead.
The default value used is `D2D_LOG_TOFILE` for the default build and `D2D_LOG_TOFILE_FULL` for a debug build.

### Change the debug level for debug builds

Before calling `decomp_2d_init`, the external code can modify the variable `decomp_debug` to change the debug level. The user can also modify this variable using the environment variable `DECOMP_2D_DEBUG`. Please note that the environment variable is processed only for debug builds. The expected value for the variable `decomp_debug` is some integer between 0 and 6, bounds included.

### List of preprocessor variables

#### DEBUG

This variable is automatically added in debug and dev builds. Extra information is printed when it is present.

#### DOUBLE_PREC

When this variable is not present, the library uses single precision. When it is present, the library uses double precision

#### SAVE_SINGLE

This variable is valid for double precision builds only. When it is present, snapshots are written in single precision.

#### PROFILER

This variable is automatically added when selecting the profiler. It activates the profiling sectionsof the code.

#### EVEN

This preprocessor variable is not valid for GPU builds. It leads to padded alltoall operations.

#### OVERWRITE

This variable leads to overwrite the input array when computing FFT. The support of this flag does not always correspond to in-place transforms, depending on the FFT backend selected, as described above.

#### HALO_DEBUG

This variable is used to debug the halo operations.

#### _GPU

This variable is automatically added in GPU builds.

#### _NCCL

This variable is valid only for GPU builds. The NVIDIA Collective Communication Library (NCCL) implements multi-GPU and multi-node communication primitives optimized for NVIDIA GPUs and Networking.

#### OCC

This variable is not supported

#### SHM

This variable is not supported.

#### SHM_DEBUG

This variable is not supported

## Optional dependencies

### FFTW

The library [fftw](http://www.fftw.org/index.html) can be used as a backend for the FFT engine. The version 3.3.10 was tested, is supported and can be downloaded [here](http://www.fftw.org/download.html). Please note that one should build fftw and decomp2d against the same compilers. For build instructions, please check [here](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html). Below is a suggestion for the compilation of the library in double precision (add `--enable-single` for a single precision build):

```
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar xzf fftw-3.3.10.tar.gz
mkdir fftw-3.3.10_tmp && cd fftw-3.3.10_tmp
../fftw-3.3.10/configure --prefix=xxxxxxx/fftw3/fftw-3.3.10_bld --enable-shared
make -j
make -j check
make install
```

### Caliper

The library [caliper](https://github.com/LLNL/Caliper) can be used to profile the execution of the code. The version 2.9.0 was tested and is supported, version 2.8.0 has also been tested and is still expected to work. Please note that one must build caliper and decomp2d against the same C/C++/Fortran compilers and MPI libray. For build instructions, please check [here](https://github.com/LLNL/Caliper#building-and-installing) and [here](https://software.llnl.gov/Caliper/CaliperBasics.html#build-and-install). Below is a suggestion for the compilation of the library using the GNU compilers:

```
git clone https://github.com/LLNL/Caliper.git caliper_github
cd caliper_github
git checkout v2.9.0
mkdir build && cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../caliper_build_2.9.0 -DWITH_FORTRAN=yes -DWITH_MPI=yes -DBUILD_TESTING=yes ../
make -j
make test
make install
```
