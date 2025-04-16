cmake_minimum_required(VERSION 3.20)                        
cmake_policy(SET CMP0074 NEW)

project(decomp2d
  LANGUAGES Fortran)
set(version 2.1.0)
enable_testing()
if (IO_BACKEND MATCHES "adios2")
  # Can be useful to also activate CXX, sometimes is needed by packages
  enable_language(C CXX)
endif (IO_BACKEND MATCHES "adios2")


set(BUILD_TARGET "mpi" CACHE STRING "Target for acceleration (mpi (default) or gpu)")
set_property(CACHE BUILD_TARGET PROPERTY STRINGS mpi gpu)

if (BUILD_TARGET MATCHES "gpu")
  option(ENABLE_OPENACC "Allow user to activate/deactivate OpenACC support" ON)
  option(ENABLE_CUDA "Allow user to activate/deactivate CUDA support" ON)
  option(ENABLE_MANAGED "Allow user to activate/deactivate automatic memory managment from NVHPC" OFF)
  option(ENABLE_NCCL "Allow user to activate/deactivate Collective Comunication NCCL" OFF)
  if (ENABLE_CUDA)
    enable_language(CUDA)
  endif()
endif(BUILD_TARGET MATCHES "gpu")

include(CMakeDependentOption)
cmake_dependent_option(ENABLE_OMP "Use OpenMP with pure MPI" ON "NOT ENABLE_CUDA" OFF)

set (ENABLE_PROFILER "OFF" CACHE STRING "Activate/select the profiler")
set_property(CACHE ENABLE_PROFILER PROPERTY STRINGS OFF caliper)
option(ENABLE_INPLACE "Enable in-place operation of FFT" OFF)
option(HALO_DEBUG "Enable halo debugging options" OFF)
option(EVEN "Padded alltoall transpose operations" OFF)
option(BUILD_SHARED_LIBS "Set ON to build a shared library" OFF)
set (IO_BACKEND "mpi" CACHE STRING "Default IO backend (mpi (default) or adios2)")
set_property(CACHE IO_BACKEND PROPERTY STRINGS mpi adios2)

set(AUTHOR "The 2DECOMP&FFT team")
set(AUTHOR_DETAILS "see https://github.com/2decomp-fft/2decomp-fft/blob/dev/CONTRIBUTORS")
set(DESCRIPTION "CMake build system for 2DECOMP&FFT")


include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/opt" CACHE PATH "..." FORCE)
endif()

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Dev Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/compilers" "${CMAKE_SOURCE_DIR}/cmake/fft")

if (BUILD_TARGET MATCHES "gpu")
  include(D2D_GPU)
endif (BUILD_TARGET MATCHES "gpu")

include(D2D_MPI)

# FFT options
if (BUILD_TARGET MATCHES "gpu")
  set(FFT_Choice "cufft" CACHE STRING "FFT for XCompact3d project (with GPU cufft is the default)")
else ()
  set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
endif (BUILD_TARGET MATCHES "gpu")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw fftw_f03 mkl cufft)

# Set now compilers flags
include(D2D_Compilers)
include(D2D_Profilers)
include(fft)

# Create a static library for the fft
add_subdirectory(src)

# Add tests
option(BUILD_TESTING "Build with tests" OFF)
if (${BUILD_TESTING})
  add_subdirectory(examples)
endif (${BUILD_TESTING})

# Add a prettify target
add_custom_target(format sh ${CMAKE_SOURCE_DIR}/scripts/format.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Write the summary 
message(STATUS "Building of  ${PROJECT_NAME} with the following options")
message( "   * BUILD TARGET   : ${BUILD_TARGET}")
if (BUILD_TARGET MATCHES "gpu")
  message( "   * NCCL           : ${ENABLE_NCCL}")
  message( "   * MANAGED MEMORY : ${ENABLE_MANAGED}")
  message( "   * CUDA ARCH      : ${CMAKE_CUDA_ARCHITECTURES}")
  message( "   * Number GPUs    : ${GPU_NUMBER}")
endif()
message( "   * BUILD EVEN     : ${EVEN}")
message( "   * PROFILER       : ${ENABLE_PROFILER}")
message( "   * OpenMP         : ${ENABLE_OMP}")
message( "   * FFT Choice     : ${FFT_Choice}")
message( "   * IO Backend     : ${IO_BACKEND}")
message( "   * IO SINGLE      : ${SINGLE_PRECISION_OUTPUT}")
