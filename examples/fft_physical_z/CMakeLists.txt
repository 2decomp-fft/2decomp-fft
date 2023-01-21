file(GLOB files_fft_c2c fft_c2c_z.f90)
file(GLOB files_fft_r2c fft_r2c_z.f90)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(fft_c2c_z ${files_fft_c2c})
add_executable(fft_r2c_z ${files_fft_r2c})

target_link_libraries(fft_c2c_z PRIVATE decomp2d)
target_link_libraries(fft_r2c_z PRIVATE decomp2d)

add_test(NAME fft_c2c_z COMMAND fft_c2c_z)
add_test(NAME fft_r2c_z COMMAND fft_r2c_z)
