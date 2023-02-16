# PROW/PCOL options
# Tests that accept PROW/PCOL as arguments will use these values.
# To simplify out of the box testing, default to 0.
set(PROW 0 CACHE STRING
  "Number of processor rows - PROWxPCOL=NP must be satisfied, 0 for autotuning")
set(PCOL 0 CACHE STRING
  "Number of processor rows - PROWxPCOL=NP must be satisfied, 0 for autotuning")
# In case decomposition is imposed force number of MPI task to be consistent
math(EXPR NUMPROCS "${PROW} * ${PCOL}")
message(STATUS "Computed NRANK: ${NUMPROCS} Max avail ${MAX_NUMPROCS}")
set(ADD_DECOMP_TO_EXAMPLE FALSE)
if (NUMPROCS GREATER "0")
  if (NUMPROCS LESS "${MAX_NUMPROCS}")
    message(STATUS "Decomposion has been imposed to ${PROW}X${PCOL}: number of MPI tasks to be used is imposed to ${NUMPROCS}")
    set(MPIEXEC_MAX_NUMPROCS "${NUMPROCS}" CACHE STRING
            "Force N to be p_row*p_col" FORCE)
  else ()
    message(STATUS "The decomposition ${PROW}x${PCOL} cannot be run. Only ${MAX_NUMPROCS} are available. Default testing is performed")
  endif()
endif ()
#string(JOIN " " TEST_ARGUMENTS "${PROW}" "${PCOL}")
set(TEST_ARGUMENTS "${PROW}" "${PCOL}")
message(STATUS "Test argument string ${TEST_ARGUMENTS}")
