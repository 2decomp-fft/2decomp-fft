# Compilers CMakeLists

set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")



if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  include(D2D_flags_gnu)
  #set(D2D_FFLAGS "-cpp -std=f2008 -ffree-line-length-none")
  #if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  #  message(STATUS "Set New Fortran basic flags")
  #  set(D2D_FFLAGS "${D2D_FFLAGS} -fallow-argument-mismatch")
  #endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  #set(D2D_FFLAGS_RELEASE "-O3 -march=native")
  #set(D2D_FFLAGS_DEBUG   "-DDEBUG -g3 -Og -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none")
  #set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -Wall -Wpedantic -Wno-unused-function -Werror -Wno-integer-division")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  message(STATUS "Setting ifort flags")
  include(D2D_flags_intel)
elseif (Fortran_COMPILER_NAME MATCHES "NAG")
  message(STATUS "Setting nagfor flags")
  include(D2D_flags_nag)
elseif (Fortran_COMPILER_NAME MATCHES "Cray")
  message(STATUS "Setting cray fortran flags")
  include(D2D_flags_cray)
elseif (Fortran_COMPILER_NAME MATCHES "NVHPC")
  message(STATUS "Setting NVHPC fortran flags")
  include(D2D_flags_nvidia)
# elseif (Fortran_COMPILER_NAME MATCHES "Flang")
#   message(STATUS "Setting Flang flags")
#   set(CMAKE_Fortran_FLAGS "-cpp -std=f2008" CACHE STRING
# 	  "Baseline FFLAGS"
# 	  FORCE)
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  include(D2D_flags_fujitsu)
else (Fortran_COMPILER_NAME MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set(D2D_FFLAGS_RELEASE "-O2")
  set(D2D_FFLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "GNU")

if (NOT FLAGS_SET)
  set(CMAKE_Fortran_FLAGS ${D2D_FFLAGS} CACHE STRING 
	"Base FFLAGS for build" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE ${D2D_FFLAGS_RELEASE} CACHE STRING
  	"Additional FFLAGS for Release (optimised) build" FORCE)
  set(CMAKE_Fortran_FLAGS_DEBUG ${D2D_FFLAGS_DEBUG} CACHE STRING
  	"Additional FFLAGS for Debug build" FORCE)
  set(CMAKE_Fortran_FLAGS_DEV ${D2D_FFLAGS_DEV} CACHE STRING
  	"Additional FFLAGS for Dev build" FORCE)
  # Add profiler
  #if (ENABLE_PROFILER)
  #  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${D2D_EXE_LINKER_FLAGS}" CACHE STRING 
  #          "Add profiler to exe" FORCE)
  #endif(ENABLE_PROFILER)

  set(FLAGS_SET 1 CACHE INTERNAL "Flags are set")
endif()

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBUG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")

if (ENABLE_INPLACE)
  add_definitions("-DOVERWRITE")
endif ()

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
