file(GLOB files_fft_c2c fft_c2c_x.f90)
file(GLOB files_fft_r2c fft_r2c_x.f90)
file(GLOB files_fft_grid fft_grid_x.f90)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(fft_c2c_x ${files_fft_c2c})
add_executable(fft_r2c_x ${files_fft_r2c})
add_executable(fft_grid_x ${files_fft_grid})

target_link_libraries(fft_c2c_x PRIVATE decomp2d)
target_link_libraries(fft_r2c_x PRIVATE decomp2d)
target_link_libraries(fft_grid_x PRIVATE decomp2d)

add_test(NAME fft_c2c_x COMMAND fft_c2c_x)
add_test(NAME fft_r2c_x COMMAND fft_r2c_x)
add_test(NAME fft_grid_x COMMAND fft_grid_x)
