cmake_minimum_required(VERSION 3.12)                        
cmake_policy(SET CMP0074 NEW)

project(2decomp LANGUAGES Fortran)
enable_testing()


set(BUILD_TARGET "mpi" CACHE STRING "Target for acceleration (mpi (default) or gpu)")
set_property(CACHE BUILD_TARGET PROPERTY STRINGS mpi gpu)

if (BUILD_TARGET MATCHES "gpu")
  option(ENABLE_OPENACC "Allow user to activate/deactivate OpenACC support" ON)
  option(ENABLE_CUDA "Allow user to activate/deactivate CUDA support" ON)
  option(ENABLE_MANAGED "Allow user to activate/deactivate automatic memory managment from NVHPC" OFF)
  option(ENABLE_NCCL "Allow user to activate/deactivate Collective Comunication NCCL" ON)

  if (ENABLE_CUDA)
    message(STATUS "Before enable CUDA")
    enable_language(CUDA)
    message(STATUS "After enable CUDA")
  endif()

endif(BUILD_TARGET MATCHES "gpu")

option(ENABLE_PROFILER "Activate the profiler" OFF)
option(ENABLE_INPLACE "Activate the profiler" OFF)
option(HALO_DEBUG "Enable halo debugging options" OFF)

set(AUTHOR "Stefano Rolfo;Charles Moulinec;Paul Bartholomew")
set(AUTHOR_DETAILS "stefano.rolfo@stfc.ac.uk;charles.moulinec@stfc.ac.uk;p.bartholomew@epcc.ed.ac.uk")
set(DESCRIPTION "Building 2decomp&fft using cmake")

message(STATUS "building ${PROJECT_NAME}")

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/mod)
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

include(D2D_Compilers)

# FFT options
if (BUILD_TARGET MATCHES "gpu")
  set(FFT_Choice "cufft" CACHE STRING "FFT for XCompact3d project (with GPU cufft is the default)")
else ()
  set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
endif (BUILD_TARGET MATCHES "gpu")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw mkl cufft)
include(fft)

# Create a static library for the fft
add_subdirectory(src)

# Add tests
add_subdirectory(examples)

# # Create an example dir with all input.i3d example files
# option(BUILD_TESTING "Build with tests" ON)
# set(test_dir "${PROJECT_BINARY_DIR}/Test")
# message(STATUS "Before test main ${test_dir}")
# if (${BUILD_TESTING})
#   file(MAKE_DIRECTORY ${test_dir})
#   include(CTest)
#   message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
#   message(STATUS "MPI EXEC: ${MPIEXEC_EXECUTABLE}")
#   file(MAKE_DIRECTORY ${test_dir}/App)
#   add_test(NAME TestAPP COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d WORKING_DIRECTORY ${test_dir}/App)
#   #if (PFUNIT_FOUND)
#   #  file(MAKE_DIRECTORY ${test_dir}/Verif_x_00)
#   #  add_test(NAME Test_dx_00 COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/verif_x_00 WORKING_DIRECTORY ${test_dir}/Verif_x_00)
#   #endif (PFUNIT_FOUND)
# endif()

