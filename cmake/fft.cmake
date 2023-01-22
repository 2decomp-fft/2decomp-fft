# FFT CMakeLists

# Look for fftw if required
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

elseif(${FFT_Choice} MATCHES "mkl")
  find_package(MKL CONFIG REQUIRED)
endif(${FFT_Choice} MATCHES "fftw")
