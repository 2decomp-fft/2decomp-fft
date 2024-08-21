# Compilers Flags for Intel
# Check is the compiler is the new ifx based on LLVM or the old ifort
if (Fortran_COMPILER_NAME MATCHES "IntelLLVM")
  set(D2D_FFLAGS "-fpp -std08 -safe-cray-ptr -g -traceback")
  set(D2D_FFLAGS_RELEASE "-O3")
  #set(D2D_FFLAGS "-fpp -std08 -xHost -heaparrays -safe-cray-ptr -g -traceback")
  #set(D2D_FFLAGS_RELEASE "-O3 -ipo")
else (Fortran_COMPILER_NAME MATCHES "IntelLLVM")
  #set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
  set(D2D_FFLAGS "-fpp -std08 -xHost -heaparrays -safe-cray-ptr -g -traceback")
  set(D2D_FFLAGS_RELEASE "-O3 -ipo")
endif (Fortran_COMPILER_NAME MATCHES "IntelLLVM")
set(D2D_FFLAGS_DEBUG   "-g -O0 -debug extended -traceback -DDEBUG")
set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -warn all,noexternal")
