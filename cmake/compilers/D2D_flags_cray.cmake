#Compilers Flags for Cray

set(D2D_FFLAGS "-eF -g -N 1023 -M878 -M296")
if (ENABLE_OPENMP_OFFLOAD)
  set(D2D_FFLAGS "${D2D_FFLAGS} -h omp")
endif()
set(D2D_FFLAGS_RELEASE "-O3")
set(D2D_FFLAGS_DEBUG   "-O0 -g")

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 19.0.0)
  message(STATUS "Setting submodules for new CrayFTN")
  set(CMAKE_Fortran_SUBMODULE_SEP ".")
  set(CMAKE_Fortran_SUBMODULE_EXT ".smod")
endif()
