cmake_minimum_required(VERSION 3.9)

project(2decomp-fft Fortran)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

find_package(MPI REQUIRED)

set(CMAKE_Fortran_FLAGS "-g -fbacktrace -fdefault-real-8 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch -cpp -O3")

file(GLOB files_decomp src/decomp_2d.f90
	               src/d2d_log.f90
                       src/glassman.f90)

#if(${FFT_Choice} MATCHES "generic")
#  file(GLOB files_fft src/fft_generic.f90)
#endif(${FFT_Choice} MATCHES "generic")

if(FFTW_FOUND)
  message (STATUS "Compiling using FFTW3")
  file(GLOB files_fft src/fft_fftw3.f90)
else(FFTW_FOUND)
  message (STATUS "Compiling using Generic FFT")
   file(GLOB files_fft src/fft_generic.f90)
endif(FFTW_FOUND)

option(DOUBLE_PRECISION "Build 2decomp_fft with double precision" ON)
if (DOUBLE_PRECISION)
  add_definitions("-DDOUBLE_PREC")
endif()

option(SINGLE_PRECISION_OUTPUT "Build 2decomp_fft with output in single precision" OFF)
if (SINGLE_PRECISION_OUTPUT)
  add_definitions("-DSAVE_SINGLE")
endif()

set(SRCFILES ${files_decomp} ${files_fft})

add_library(2decomp-fft STATIC ${SRCFILES})
if (MPI_FOUND)
  target_link_libraries(2decomp-fft PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(FFTW_FOUND)
  target_include_directories(2decomp-fft PRIVATE ${FFTW_INCLUDE_DIRS})
  #target_include_directories(2decomp-fft PRIVATE ${PKG_FFTW_LIBRARY_DIRS})
endif(FFTW_FOUND)

#install(TARGETS decomp2d
#        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
#)

#add_subdirectory(examples)
#include_directories(examples)
