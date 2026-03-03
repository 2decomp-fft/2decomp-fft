file(GLOB files_test interp.f90)
set(case "interp")
set(case_dir "${test_dir}/${case}")
CreateMPITest("${case_dir}" "${case}" "${files_test}" "" "" "" "1")
