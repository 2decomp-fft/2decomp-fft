#
# 2decomp CMake configuration.
# This is used by other packages to configure themselves against 2decomp&fft.
#

# Compute installation prefix relative to this file
get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_prefix "${_dir}/../.." ABSOLUTE)

# Import the targets
if (EXISTS ${_prefix}/lib )
  message(STATUS "Found decomp2d under lib")
  include("${_prefix}/lib/decomp2d-targets.cmake")
else()
  message(STATUS "Not Found decomp2d under lib, we'll use lib64")
  include("${_prefix}/lib64/decomp2d-targets.cmake")
endif()

set(2decomp_INCLUDE_DIR "${_prefix}/include")

