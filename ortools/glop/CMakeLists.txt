file(GLOB_RECURSE _SRCS "*.h" "*.cc")
set(NAME ${PROJECT_NAME}_glop)

# Will be merge in libortools.so
add_library(${NAME} STATIC ${_SRCS})
#add_library(${NAME} OBJECT ${_SRCS})
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_STANDARD 11)
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
set_target_properties(${NAME} PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
set_target_properties(${NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(${NAME} PRIVATE
	${PROJECT_SOURCE_DIR}
	${PROJECT_BINARY_DIR})
target_link_libraries(${NAME} PRIVATE Protobuf gflags glog Cbc ${PROJECT_NAME}::proto)
add_library(${PROJECT_NAME}::glop ALIAS ${NAME})
#add_dependencies(${NAME} ${PROJECT_NAME}::proto)
