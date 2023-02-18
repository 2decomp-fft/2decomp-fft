# MPI CMakeLists

find_package(MPI REQUIRED)
# Stop if there is no MPI_Fortran_Compiler
if (MPI_Fortran_COMPILER)
    message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
    message(STATUS "MPI_VERSION found: ${MPI_VERSION}")
    # Try to guess the MPI type to adapt compilation flags if necessary
    string(FIND "${MPI_Fortran_COMPILER}" "mpich" pos)
    if(pos GREATER_EQUAL "0")
      set(FIND_MPICH TRUE)
      message(STATUS "MPI is MPICH type")
    endif()
    string(FIND "${MPI_Fortran_COMPILER}" "openmpi" pos)
    if(pos GREATER_EQUAL "0")
      set(FIND_OMPI TRUE)
      message(STATUS "MPI is openMPI type")
    endif()
else (MPI_Fortran_COMPILER)
    message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)
# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
    message(STATUS "MPI FOUND: ${MPI_FOUND}")
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
    if (NOT MPI_NUMPROCS_SET)
      message(STATUS "Reset the number of ranks to 1")
      set(MPIEXEC_MAX_NUMPROCS "1" CACHE STRING
          "Set the initial value to 1 rank" FORCE)
      set(MPI_NUMPROCS_SET 1 CACHE INTERNAL "MPI Ranks set")
      # Force the mpirun to be coherent with the mpifortran
      string(REGEX REPLACE "mpif90" "mpirun" PATH_TO_MPIRUN "${MPI_Fortran_COMPILER}")
      string(REPLACE "mpiifort" "mpirun" PATH_TO_MPIRUN "${PATH_TO_MPIRUN}")
      message(STATUS "Path to mpirun ${PATH_TO_MPIRUN}")
      set(MPIEXEC_EXECUTABLE "${PATH_TO_MPIRUN}" CACHE STRING
          "Force MPIRUN to be consistent with MPI_Fortran_COMPILER" FORCE)
    endif()
else (MPI_FOUND)
    message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
    message(STATUS "Code will compile but performaces can be compromised")
    message(STATUS "Using a CMake vers > 3.10 should solve the problem")
    message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)

# PROW/PCOL options
# Tests that accept PROW/PCOL as arguments will use these values.
# To simplify out of the box testing, default to 0.
set(PROW 0 CACHE STRING
  "Number of processor rows - PROWxPCOL=NP must be satisfied, 0 for autotuning")
set(PCOL 0 CACHE STRING
  "Number of processor rows - PROWxPCOL=NP must be satisfied, 0 for autotuning")
