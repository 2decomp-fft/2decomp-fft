# FFT CMakeLists

# Look for fftw if required
# The download findFFTW code is based on the README at
# https://github.com/egpbos/findFFTW
#
# Copyright (c) 2015, Wenzel Jakob; 2017, Patrick Bos
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
message(STATUS "SET UP FFT")
if(${FFT_Choice} MATCHES "fftw")
  configure_file(${CMAKE_SOURCE_DIR}/cmake/fft/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
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
elseif(${FFT_Choice} MATCHES "cufft")
  message(STATUS "Enable cuFFT")
  if (ENABLE_CUDA)
    set(CUFFT_FOUND TRUE)
  endif()
endif(${FFT_Choice} MATCHES "fftw")
