file(GLOB files_test init_test.f90)
set(case "init_test")
set(case_dir "${test_dir}/${case}")
CreateMPITest("${case_dir}" "${case}" "${files_test}" "" "" "" "1")
