cmake_minimum_required(VERSION 3.3)
project(kdv5-example)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")

add_subdirectory(kdv5 EXCLUDE_FROM_ALL)

add_executable(kdv5-example-1 kdv5_example_1.cpp)
target_link_libraries(kdv5-example-1 kdv5)

add_custom_command(TARGET kdv5-example-1 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        $<TARGET_FILE_DIR:lapacke>
        $<TARGET_FILE_DIR:kdv5-example-1>)     
		
add_executable(kdv5-example-2 kdv5_example_2.cpp)
target_link_libraries(kdv5-example-2 kdv5)

add_custom_command(TARGET kdv5-example-2 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        $<TARGET_FILE_DIR:lapacke>
        $<TARGET_FILE_DIR:kdv5-example-2>)   		