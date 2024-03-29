# Flags for GNU compiler
set(D2D_FFLAGS "-cpp -std=f2008 -ffree-line-length-none")
if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  message(STATUS "Set New Fortran basic flags")
  set(D2D_FFLAGS "${D2D_FFLAGS} -fallow-argument-mismatch")
  set(D2D_GNU10 TRUE)
else (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(D2D_GNU10 FALSE)
endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
set(D2D_FFLAGS_RELEASE "-O3 -march=native")
set(D2D_FFLAGS_DEBUG   "-DDEBUG -g3 -Og -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none")
if (((${D2D_MPI_FAMILY} STREQUAL "MPICH") OR (${D2D_MPI_FAMILY} STREQUAL "Unknown")) AND D2D_GNU10)
  set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -Wall -Wno-unused-function -Wno-integer-division")
else()
  set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG} -Wall -Wpedantic -Wno-unused-function -Werror -Wno-integer-division")
endif()
if ((${D2D_MPI_FAMILY} STREQUAL "OMPI") AND (${FFT_Choice} MATCHES "generic"))
  set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEV} -Wimplicit-procedure -Wimplicit-interface")
endif()
