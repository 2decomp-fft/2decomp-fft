file(GLOB files_decomp factor.f90
	               decomp_2d_constants.f90
		       decomp_2d_mpi.f90
	               decomp_2d.f90
		       log.f90
		       io.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

if(${BUILD_TARGET} MATCHES "gpu")
  list(APPEND files_decomp decomp_2d_cumpi.f90)
  if(ENABLE_NCCL)
    list(APPEND files_decomp decomp_2d_nccl.f90)
  endif()
endif(${BUILD_TARGET} MATCHES "gpu")

if(FFTW_FOUND)
  if (${FFT_Choice} MATCHES "fftw_f03")
    message (STATUS "Compiling using FFTW3 F2003 interface")
    file(GLOB files_fft fft_fftw3_f03.f90)
  else()
    message (STATUS "Compiling using FFTW3 legacy interface")
    file(GLOB files_fft fft_fftw3.f90)
  endif()
elseif(MKL_FOUND)
  message (STATUS "Compiling using MKL")
  file(GLOB files_fft $ENV{MKLROOT}/include/mkl_dfti.f90 fft_mkl.f90)
elseif(CUFFT_FOUND)
  message (STATUS "Compiling using cuFFT")
  file(GLOB files_fft fft_cufft.f90)
else(FFTW_FOUND)
  message (STATUS "Compiling using Generic FFT")
  file(GLOB files_fft glassman.f90 fft_generic.f90)
endif(FFTW_FOUND)
list(APPEND files_fft fft_log.f90)

if (ENABLE_PROFILER)
  file(GLOB prof_files profiler_${ENABLE_PROFILER}.f90)
endif()

set(SRCFILES ${files_decomp} ${files_fft} ${prof_files})

add_library(decomp2d ${SRCFILES})
target_include_directories(decomp2d PUBLIC
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if (MPI_FOUND)
  target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(FFTW_FOUND)
  message(STATUS "FFTW_DOUBLE     : ${FFTW_DOUBLE_LIB}")
  message(STATUS "FFTW_SINGLE     : ${FFTW_FLOAT_LIB}")
  target_include_directories(decomp2d PRIVATE ${FFTW_INCLUDE_DIRS})
elseif(MKL_FOUND)
  target_compile_options(decomp2d PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_include_directories(decomp2d PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
endif(FFTW_FOUND)

if (ENABLE_PROFILER)
  target_compile_definitions(decomp2d PUBLIC -DPROFILER)

  if (caliper_FOUND)
    target_include_directories(decomp2d PUBLIC ${caliper_INCLUDE_DIR}/caliper/fortran)
    target_link_libraries(decomp2d PRIVATE caliper)
  endif()
endif()
if (HALO_DEBUG)
  target_compile_definitions(decomp2d PUBLIC -DHALO_DEBUG)
endif()

install(TARGETS decomp2d
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.mod")
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.smod")
      
