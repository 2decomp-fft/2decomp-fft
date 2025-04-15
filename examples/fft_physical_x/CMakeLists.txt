file(GLOB files_fft_c2c fft_c2c_x.f90)
set(case "fft_c2c_x")
CreateMPITest(${test_dir} ${case} ${files_fft_c2c} "1")

file(GLOB files_fft_r2c fft_r2c_x.f90)
set(case "fft_r2c_x")
CreateMPITest(${test_dir} ${case} ${files_fft_r2c} "1")

file(GLOB files_fft_grid fft_grid_x.f90)
set(case "fft_grid_x")
CreateMPITest(${test_dir} ${case} ${files_fft_grid} "1")

file(GLOB files_fft_c2c_skip fft_c2c_x_skip.f90)
set(case "fft_c2c_skip_x")
CreateMPITest(${test_dir} ${case} ${files_fft_c2c_skip} "1")

file(GLOB files_fft_r2c_skip fft_r2c_x_skip.f90)
set(case "fft_r2c_skip_x")
CreateMPITest(${test_dir} ${case} ${files_fft_r2c_skip} "1")

file(GLOB files_fft_r2c_skip_warn fft_r2c_x_skip_errorMsg.f90)
set(case "fft_r2c_skip_x_warn")
CreateMPITest(${test_dir} ${case} ${files_fft_r2c_skip_warn} "1")



