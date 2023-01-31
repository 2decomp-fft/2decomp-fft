file(GLOB files_test halo_test.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(halo_test ${files_test})
target_link_libraries(halo_test PRIVATE decomp2d)
add_test(NAME halo_test COMMAND ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:halo_test>)
