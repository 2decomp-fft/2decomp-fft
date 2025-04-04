# Compilers CMakeLists

set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")



if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  include(D2D_flags_gnu)
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
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  include(D2D_flags_fujitsu)
elseif (Fortran_COMPILER_NAME MATCHES "Flang")
  message(STATUS "Setting Flang flags")
  include(D2D_flags_flang)
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

# Padded MPI alltoall transpose operations (invalid for GPU)
if (EVEN)
	#if (BUILD_TARGET MATCHES "gpu")
  	#  message(FATAL_ERROR "The GPU build is not compatible with padded alltoall")
  	#endif ()
  add_definitions("-DEVEN")
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

if (IO_BACKEND MATCHES "mpi")
  message(STATUS "Using mpi (default) IO backend")
elseif (IO_BACKEND MATCHES "adios2")
  message(STATUS "Using ADIOS2 IO backend")
  find_package(adios2 REQUIRED)
  if (NOT ADIOS2_HAVE_MPI)
    message(FATAL_ERROR "MPI support is missing in the provided ADIOS2 build")
  endif (NOT ADIOS2_HAVE_MPI)
  if (NOT ADIOS2_HAVE_Fortran)
    message(FATAL_ERROR "Fortran support is missing in the provided ADIOS2 build")
  endif (NOT ADIOS2_HAVE_Fortran)
else (IO_BACKEND MATCHES "mpi")
  message(FATAL_ERROR "Invalid value for CMake variable IO_BACKEND")
endif (IO_BACKEND MATCHES "mpi")
