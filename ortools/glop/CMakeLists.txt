file(GLOB _SRCS "*.h" "*.cc")
set(NAME ${PROJECT_NAME}_glop)

# Will be merge in libortools.so
#add_library(${NAME} STATIC ${_SRCS})
add_library(${NAME} OBJECT ${_SRCS})
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_STANDARD 11)
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
set_target_properties(${NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(${NAME} PRIVATE
	${PROJECT_SOURCE_DIR}
	${PROJECT_BINARY_DIR}
	$<TARGET_PROPERTY:Protobuf,INTERFACE_INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:gflags,INTERFACE_INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:glog,INTERFACE_INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:Cbc,INTERFACE_INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:${PROJECT_NAME}::proto,INTERFACE_INCLUDE_DIRECTORIES>)
#target_link_libraries(${NAME} PRIVATE Protobuf gflags glog Cbc ${PROJECT_NAME}::proto)
add_dependencies(${NAME} Protobuf gflags glog Cbc ${PROJECT_NAME}::proto)
add_library(${PROJECT_NAME}::glop ALIAS ${NAME})
