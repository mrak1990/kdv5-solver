cmake_minimum_required(VERSION 3.3)
project(kdv5-solver CXX Fortran)

if (NOT CMAKE_BUILD_TYPE) 
      set(CMAKE_BUILD_TYPE Release) 
endif() 

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")

add_subdirectory(source)
add_subdirectory(lib)
add_subdirectory(tests)