# Compilers Flags for Intel

set(D2D_FFLAGS "-fpp -std08 -xHost -heaparrays -safe-cray-ptr -g -traceback")
set(D2D_FFLAGS_RELEASE "-O3 -ipo")
set(D2D_FFLAGS_DEBUG   "-g -O0 -debug extended -traceback -DDEBUG")
set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -warn all,noexternal")
#set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
