file(GLOB files_dtt_x fftw_f03_dtt_x.f90)
file(GLOB files_dtt_z fftw_f03_dtt_z.f90)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${FFTW_INCLUDE_DIRS})

add_executable(dtt_x ${files_dtt_x})
add_executable(dtt_z ${files_dtt_z})

target_link_libraries(dtt_x PRIVATE decomp2d examples_utils)
target_link_libraries(dtt_z PRIVATE decomp2d examples_utils)

# Run the test(s)
set(run_dir "${test_dir}/dtt_x")
message(STATUS "Example dir ${run_dir}")
file(MAKE_DIRECTORY ${run_dir})
add_test(NAME dtt_x COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:dtt_x> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
set(run_dir "${test_dir}/dtt_z")
message(STATUS "Example dir ${run_dir}")
file(MAKE_DIRECTORY ${run_dir})
add_test(NAME dtt_z COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:dtt_z> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
