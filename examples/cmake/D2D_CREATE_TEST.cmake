macro(CreateMPITest run_dir
                    case
		    app_src
		    exe_dep
		    defs
		    files
                    OMP_THREADS)
  message(STATUS "Add Verification Test (MPI run) ${case}")
  # Create working directory
  file(MAKE_DIRECTORY ${run_dir})
  # Copy additional files
  set(loc_files "")
  list(APPEND loc_files ${files})
  if (BUILD_TARGET MATCHES "gpu")
    list(APPEND loc_files "bind.sh")
  endif()
  foreach(ff IN LISTS loc_files)
    file(COPY ${ff} DESTINATION ${run_dir})
  endforeach()
  # Add includes directories and create executable
  include_directories(${CMAKE_SOURCE_DIR}/src)
  if (FFT_Choice MATCHES "fftw_f03")
    include_directories(${FFTW_INCLUDE_DIRS})
  endif()
  add_executable(${case} ${app_src})
  # Add dependecies
  set(loc_exe_dep "")
  list(APPEND loc_exe_dep ${exe_dep})
  foreach(dep IN LISTS loc_exe_dep)
    add_dependencies(${case} ${exe_dep})
  endforeach()
  # Add definitions
  set(loc_defs "")
  list(APPEND loc_defs ${defs})
  foreach(def IN LISTS loc_defs)
    add_compile_definitions(def)
  endforeach()
  # Linking
  target_link_libraries(${case} PRIVATE decomp2d examples_utils)
  if (OPENMP_FOUND AND ENABLE_OMP)
    target_link_libraries(${case} PRIVATE OpenMP::OpenMP_Fortran)
    if (Fortran_COMPILER_NAME MATCHES "Cray")
      target_link_options(${case} PRIVATE -h omp)
    endif()
  endif()
  # Install
  install(TARGETS ${case} DESTINATION ${run_dir})
  # Run the test
  if (BUILD_TARGET MATCHES "gpu")
    add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROCS} ./bind.sh $<TARGET_FILE:${case}> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
  else ()  
    add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROCS} $<TARGET_FILE:${case}> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
    if (OPENMP_FOUND AND ENABLE_OMP)
      set_tests_properties(${case} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=${OMP_THREADS}")
    endif (OPENMP_FOUND AND ENABLE_OMP)
  endif()
endmacro()
