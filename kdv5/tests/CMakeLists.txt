add_executable(kdv5-test main.cpp l2_norm.cpp approximation_order.cpp gtest_ext.h)
target_link_libraries(kdv5-test kdv5 eigen service gtest)
target_include_directories(kdv5-test PUBLIC ${CMAKE_SOURCE_DIR}/include)