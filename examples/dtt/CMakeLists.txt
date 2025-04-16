file(GLOB files_dtt_x fftw_f03_dtt_x.f90)
set(case "dtt_x")
set(case_dir "${test_dir}/${case}")
CreateMPITest("${case_dir}" "${case}" "${files_dtt_x}" "" "" "" "1")

file(GLOB files_dtt_z fftw_f03_dtt_z.f90)
set(case "dtt_z")
set(case_dir "${test_dir}/${case}")
CreateMPITest("${case_dir}" "${case}" "${files_dtt_z}" "" "" "" "1")

#include_directories(${FFTW_INCLUDE_DIRS})

