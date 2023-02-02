# Compilers CMakeLists

set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")

if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  set(D2D_FFLAGS "-cpp -std=f2008 -ffree-line-length-none")
  if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
    message(STATUS "Set New Fortran basic flags")
    set(D2D_FFLAGS "${D2D_FFLAGS} -fallow-argument-mismatch")
  endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(D2D_FFLAGS_RELEASE "-O3 -march=native")
  set(D2D_FFLAGS_DEBUG   "-DDEBUG -g3 -Og -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none")
  set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -Wall -Wpedantic -Wno-unused-function -Werror -Wno-integer-division")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  message(STATUS "Setting ifort flags")
  set(D2D_FFLAGS "-fpp -std08 -xHost -heaparrays -safe-cray-ptr -g -traceback")
  set(D2D_FFLAGS_RELEASE "-O3 -ipo")
  set(D2D_FFLAGS_DEBUG   "-g -O0 -debug extended -traceback -DDEBUG")
  set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -warn all,noexternal")
  #set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
elseif (Fortran_COMPILER_NAME MATCHES "NAG")
  message(STATUS "Setting nagfor flags")
  set(D2D_FFLAGS "-fpp")
  set(D2D_FFLAGS_RELEASE "-O3")
  set(D2D_FFLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "Cray")
  message(STATUS "Setting cray fortran flags")
  set(D2D_FFLAGS "-eF -g -N 1023")
  set(D2D_FFLAGS_RELEASE "-O3")
  set(D2D_FFLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "PGI")
  message(STATUS "Setting PGI fortran flags")
  if (ENABLE_OPENACC) 
    if (OPENACC_TARGET MATCHES "gpu")
      set(D2D_FFLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -acc -target=gpu")
    else(OPENACC_TARGET MATCHES "gpu")
      set(D2D_FFLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -acc -target=multicore")
    endif(OPENACC_TARGET MATCHES "gpu")
  else()
    set(D2D_FFLAGS "-cpp -Mfree -Kieee -Minfo=accel -g")
  endif()
  #set(CMAKE_Fortran_FLAGS "-fast -cpp -Mfree -Kieee -Minfo=accel -g ")
  set(D2D_FFLAGS_RELEASE "-fast -O3")
  set(D2D_FFLAGS_DEBUG   "-O0")
elseif (Fortran_COMPILER_NAME MATCHES "NVHPC")
  message(STATUS "Setting NVHPC fortran flags")
  if (ENABLE_OPENACC) 
    if (OPENACC_TARGET MATCHES "gpu")
      #set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -stdpar=gpu -acc -target=gpu -Minstrument")
      set(D2D_FFLAGS "${CMAKE_Fortran_FLAGS} -cpp -Mfree -Kieee -Minfo=accel -stdpar=gpu -gpu=cc${CUDA_ARCH_FIRST},managed -acc -target=gpu")
      #add_definitions("-lnvhpcwrapnvtx")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lnvhpcwrapnvtx")
    else()
      set(D2D_FFLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -stdpar -acc -target=multicore")
    endif(OPENACC_TARGET MATCHES "gpu")
  else()
    set(D2D_FFLAGS "-cpp -Mfree -Kieee -Minfo=accel -g")
  endif()
  set(D2D_FFLAGS_RELEASE "-fast -O3")
  set(D2D_FFLAGS_DEBUG   "-O0 -g -DDEBUG")
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  set(D2D_FFLAGS "-Cpp")
  set(D2D_FFLAGS_RELEASE "-O3")
  set(D2D_FFLAGS_DEBUG   "-O0")
# elseif (Fortran_COMPILER_NAME MATCHES "Flang")
#   message(STATUS "Setting Flang flags")
#   set(CMAKE_Fortran_FLAGS "-cpp -std=f2008" CACHE STRING
# 	  "Baseline FFLAGS"
# 	  FORCE)
else (Fortran_COMPILER_NAME MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set(D2D_FFLAGS_RELEASE "-O2")
  set(D2D_FFLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "GNU")

set(CMAKE_Fortran_FLAGS ${D2D_FFLAGS} CACHE STRING
  "Baseline compiler FFLAGS"
  FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE ${D2D_FFLAGS_RELEASE} CACHE STRING
  "Additional FFLAGS for Release (optimised) build"
  FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG ${D2D_FFLAGS_DEBUG} CACHE STRING
  "Additional FFLAGS for Debug build"
  FORCE)
set(CMAKE_Fortran_FLAGS_DEV ${D2D_FFLAGS_DEV} CACHE STRING
  "Additional FFLAGS for Dev build"
  FORCE)

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBUG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")

execute_process(
  COMMAND git describe --tag --long --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DVERSION=\"${GIT_VERSION}\"")
option(DOUBLE_PRECISION "Build Xcompact with double precision" ON)
if (DOUBLE_PRECISION)
  add_definitions("-DDOUBLE_PREC")
endif()

option(SINGLE_PRECISION_OUTPUT "Build XCompact with output in single precision" OFF)
if (SINGLE_PRECISION_OUTPUT)
  add_definitions("-DSAVE_SINGLE")
endif()
