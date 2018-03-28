add_library(kdv5 STATIC kdv_5.cpp dgbsv.cpp exp_parameters.cpp model_parameters.cpp coefficients.cpp mesh_props.cpp)
target_include_directories(kdv5 PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(kdv5 lapacke eigen)

add_library(service STATIC service.cpp)
target_include_directories(service PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(service eigen)