file(GLOB files_test init_test.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(init_test ${files_test})
include_directories(${CMAKE_SOURCE_DIR}/src)
target_link_libraries(init_test PRIVATE decomp2d)
add_test(NAME init_test COMMAND init_test)
