function(closest_power_of_2 number result )
  set(value ${number}) 
  set(power 1) 
  while(${power} LESS_EQUAL ${value})
    math(EXPR power "${power} * 2")
  endwhile()
  math(EXPR power "${power} / 2")
  set(${result} ${power} PARENT_SCOPE)                                                                          
endfunction()

# MPI CMakeLists
find_package(MPI REQUIRED COMPONENTS Fortran)
set(D2D_MPI_FAMILY "Unknown")

# adios2 IO backend requires C and C++ MPI components
if (IO_BACKEND MATCHES "adios2")
    find_package(MPI REQUIRED COMPONENTS C CXX)
endif (IO_BACKEND MATCHES "adios2")

# Stop if there is no MPI_Fortran_Compiler
if (MPI_Fortran_COMPILER)
    message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
    message(STATUS "MPI_VERSION found: ${MPI_VERSION}")
    # Try to guess the MPI type to adapt compilation flags if necessary
    string(FIND "${MPI_Fortran_COMPILER}" "mpich" pos)
    if(pos GREATER_EQUAL "0")
      set(D2D_MPI_FAMILY "MPICH")
      message(STATUS "MPI is MPICH type")
    endif()
    string(FIND "${MPI_Fortran_COMPILER}" "openmpi" pos)
    if(pos GREATER_EQUAL "0")
      set(D2D_MPI_FAMILY "OMPI")
      message(STATUS "MPI is openMPI type")
    endif()

    if (${D2D_MPI_FAMILY} STREQUAL "Unknown")
      execute_process(COMMAND ${MPI_Fortran_COMPILER} "-show"
	OUTPUT_VARIABLE mpi_show
	ERROR_QUIET)
      
      string(FIND "${mpi_show}" "openmpi" pos)
      if(pos GREATER_EQUAL "0")
	set(D2D_MPI_FAMILY "OMPI")
      endif()

      string(FIND "${mpi_show}" "mpich" pos)
      if(pos GREATER_EQUAL "0")
	set(D2D_MPI_FAMILY "MPICH")
      endif()
    endif()

    message(STATUS "MPI Compiler family: ${D2D_MPI_FAMILY}")
else (MPI_Fortran_COMPILER)
    message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)
# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
  message(STATUS "MPI FOUND: ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
  if (NOT MPI_NUMPROCS_SET)
    # Save the Maximim number of MPI ranks on the system
    if (ENABLE_CUDA) 
       set(NP ${GPU_NUMBER})
    else ()
       set(NP ${MPIEXEC_MAX_NUMPROCS})
    endif()
    # For even we'll test with a power of 2 number of MPI RANKS
    if (EVEN)
      closest_power_of_2(${NP} NP)
    endif()
    message(STATUS "NUMBER OF PROCS USED FOR TESTING ${NP}")
    set(MPI_NUMPROCS ${NP} CACHE STRING "SAVE NRANKS FOR MPIRUN" FORCE)
    set(MPI_NUMPROCS_SET 1 CACHE INTERNAL "MPI Ranks set" FORCE)
    # Force the mpirun to be coherent with the mpifortran
    string(REGEX REPLACE "mpif90" "mpirun" PATH_TO_MPIRUN "${MPI_Fortran_COMPILER}")
    string(REPLACE "mpiifort" "mpirun" PATH_TO_MPIRUN "${PATH_TO_MPIRUN}")
    string(REPLACE "mpiifx" "mpirun" PATH_TO_MPIRUN "${PATH_TO_MPIRUN}")
    message(STATUS "Path to mpirun ${PATH_TO_MPIRUN}")
    set(MPIEXEC_EXECUTABLE "${PATH_TO_MPIRUN}" CACHE STRING
        "Force MPIRUN to be consistent with MPI_Fortran_COMPILER" FORCE)
  endif()
else (MPI_FOUND)
  message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
  message(STATUS "Code will compile but performaces can be compromised")
  message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)
