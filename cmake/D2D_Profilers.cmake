# Profilers CMakeLists

if (ENABLE_PROFILER)

  if (ENABLE_PROFILER MATCHES "caliper")
    enable_language(CXX)
    find_package(caliper REQUIRED)
  endif()
  
endif()
