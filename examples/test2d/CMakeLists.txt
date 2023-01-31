file(GLOB files_test test2d.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(test2d ${files_test})
target_link_libraries(test2d PRIVATE decomp2d)
add_test(NAME test2d COMMAND ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:test2d>)
