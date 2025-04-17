file(GLOB files_test halo_test.f90)
set(case "halo_test")
set(case_dir "${test_dir}/${case}")
CreateMPITest("${case_dir}" "${case}" "${files_test}" "" "" "" "1")
