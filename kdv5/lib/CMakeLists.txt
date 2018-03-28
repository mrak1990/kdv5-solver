# lapack/lapacke
set(LAPACKE ON CACHE BOOL "Build LAPACKE" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
add_subdirectory(lapack-3.8.0)

# eigen
add_library(eigen INTERFACE IMPORTED GLOBAL)
set_target_properties(eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/lib/eigen-eigen-5a0156e40feb)

# gtest
set(BUILD_GMOCK ON CACHE BOOL "Builds the googlemock subproject" FORCE)
set(BUILD_GTEST OFF CACHE BOOL "Builds the googletest subproject" FORCE)
add_subdirectory(googletest-master)