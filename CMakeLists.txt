cmake_minimum_required(VERSION 3.12)                        
cmake_policy(SET CMP0074 NEW)

project(2decomp LANGUAGES Fortran)
enable_testing()

option(ENABLE_OPENACC "Allow user to activate/deactivate OpenACC support" OFF)
option(ENABLE_CUDA "Allow user to activate/deactivate CUDA support" OFF)

if (ENABLE_CUDA)
  enable_language(CUDA)
endif()

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

# Find the modules included with Xcompact
#SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Dev Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if (ENABLE_OPENACC)
  include(FindOpenACC)
  if(OpenACC_Fortran_FOUND)
    message(STATUS "OpenACC for Fotran Compiler Found, version ${OpenACC_Fortran_VERSION_MAJOR}.${OpenACC_Fortran_VERSION_MINOR}")
  else()
    message(ERROR_CRITICAL "No OpenACC support detected")
  endif()
  set(OPENACC_TARGET "gpu" CACHE STRING "Target for acceleration (gpu (default) or multicore)")
  set_property(CACHE OPENACC_TARGET PROPERTY STRINGS gpu multicore)
endif()

if (OPENACC_TARGET MATCHES "gpu")
  find_package(CUDAToolkit REQUIRED)
  if(${CMAKE_VERSION} VERSION_LESS_EQUAL "3.13.4")
    cuda_select_nvcc_arch_flags(ARCH_FLAGS "Auto") # optional argument for arch to add
    message(STATUS "ARCH_FLAGS = ${ARCH_FLAGS}")
    string(REPLACE "-gencode;" "--generate-code=" ARCH_FLAGS "${ARCH_FLAGS}")
    string(APPEND CMAKE_CUDA_FLAGS "${ARCH_FLAGS}")
    message(STATUS "ARCH_FLAGS WITH CUDA = ${ARCH_FLAGS}")
  else()
    include(FindCUDA/select_compute_arch)
    CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
    string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
    string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
    string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
    SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
    set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
    message(STATUS "CUDA_ARCHITECTURES ${CUDA_ARCH_LIST}")
    list(GET CUDA_ARCH_LIST 0 CUDA_ARCH_FIRST)
    message(STATUS "CUDA_ARCH_FIRST ${CUDA_ARCH_FIRST}")
  endif()
  message(STATUS "CUDA_LIBRARIES ${CUDA_LIBRARIES}")
endif()

set(NVIDIA_PROFILER "Use the NVidia profiles" ON)
set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")
if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2008 -ffree-line-length-none")
  if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
    message(STATUS "Set New Fortran basic flags")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-DDEBUG -g3 -Og -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_DEV     "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -Wpedantic -Wno-unused-function -Werror -Wno-integer-division")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  message(STATUS "Setting ifort flags")
  set(CMAKE_Fortran_FLAGS "-fpp -std08 -xHost -heaparrays -safe-cray-ptr -g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -debug extended -traceback -DDEBUG")
  set(CMAKE_Fortran_FLAGS_DEV     "${CMAKE_Fortran_FLAGS_DEBUG} -warn all,noexternal")
  #set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
elseif (Fortran_COMPILER_NAME MATCHES "NAG")
  message(STATUS "Setting nagfor flags")
  set(CMAKE_Fortran_FLAGS "-fpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "Cray")
  message(STATUS "Setting cray fortran flags")
  set(CMAKE_Fortran_FLAGS "-eF -g -N 1023")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "PGI")
  message(STATUS "Setting PGI fortran flags")
  if (ENABLE_OPENACC) 
    if (OPENACC_TARGET MATCHES "gpu")
      set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -acc -target=gpu")
    else(OPENACC_TARGET MATCHES "gpu")
      set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -acc -target=multicore")
    endif(OPENACC_TARGET MATCHES "gpu")
  else()
    set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g")
  endif()
  #set(CMAKE_Fortran_FLAGS "-fast -cpp -Mfree -Kieee -Minfo=accel -g ")
  set (CMAKE_Fortran_FLAGS_RELEASE "-fast -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0")
elseif (Fortran_COMPILER_NAME MATCHES "NVHPC")
  message(STATUS "Setting NVHPC fortran flags")
  if (ENABLE_OPENACC) 
    if (OPENACC_TARGET MATCHES "gpu")
      #set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -stdpar=gpu -acc -target=gpu -Minstrument")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -Mfree -Kieee -Minfo=accel -stdpar=gpu -gpu=cc${CUDA_ARCH_FIRST},managed -acc -target=gpu")
      #add_definitions("-lnvhpcwrapnvtx")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lnvhpcwrapnvtx")
    else()
      set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g -stdpar -acc -target=multicore")
    endif(OPENACC_TARGET MATCHES "gpu")
  else()
    set(CMAKE_Fortran_FLAGS "-cpp -Mfree -Kieee -Minfo=accel -g")
  endif()
  set (CMAKE_Fortran_FLAGS_RELEASE "-fast -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -DDEBUG")
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  set (CMAKE_Fortran_FLAGS "-Cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
else (Fortran_COMPILER_NAME MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "GNU")

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBUG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")


find_package(MPI REQUIRED)
# Stop if there is no MPI_Fortran_Compiler
if (MPI_Fortran_COMPILER)
    message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
else (MPI_Fortran_COMPILER)
    message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)
# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
    message(STATUS "MPI FOUND: ${MPI_FOUND}")
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
else (MPI_FOUND)
    message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
    message(STATUS "Code will compile but performaces can be compromised")
    message(STATUS "Using a CMake vers > 3.10 should solve the problem")
    message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)

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

# FFT options
set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw mkl)

# Look for fftw is required
if(${FFT_Choice} MATCHES "fftw")
  configure_file(cmake/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
      else()
      message("CMake step for findFFTW completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
  endif()
  
  set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
  
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
  find_package(FFTW)
  message(STATUS "FFTW_FOUND     : ${FFTW_FOUND}")
  message(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
  message(STATUS "FFTW_INCLUDE   : ${FFTW_INCLUDE_DIRS}")

  #add_definitions("-lfftw3 -lfftw3f")

endif(${FFT_Choice} MATCHES "fftw")

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

