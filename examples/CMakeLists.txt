# Set the decomposition and resolution
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/examples/cmake")
include(D2D_Set_Decomposition)
include(D2D_Set_Resolution)
# Create a folder directory for testing
set(test_dir "${PROJECT_BINARY_DIR}/RunTests")
file(MAKE_DIRECTORY ${test_dir})
# add subdirectories
add_subdirectory(init_test)
add_subdirectory(test2d)
add_subdirectory(fft_physical_x)
add_subdirectory(fft_physical_z)
add_subdirectory(halo_test)
add_subdirectory(io_test)

# Set real/complex tests
set(COMPLEX_TESTS "OFF" CACHE STRING "Enables complex numbers for tests that support it")
