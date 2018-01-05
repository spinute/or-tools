if (NOT BUILD_CXX)
	return()
endif()

include(utils)
set_version(VERSION)
project(ortools LANGUAGES CXX VERSION ${VERSION})
message(STATUS "ortools version: ${PROJECT_VERSION}")

# config options
if (MSVC)
	# /wd4005  macro-redefinition
	# /wd4068  unknown pragma
	# /wd4244  conversion from 'type1' to 'type2'
	# /wd4267  conversion from 'size_t' to 'type2'
	# /wd4800  force value to bool 'true' or 'false' (performance warning)
	add_compile_options(/W3 /WX /wd4005 /wd4068 /wd4244 /wd4267 /wd4800)
	add_definitions(/DNOMINMAX /DWIN32_LEAN_AND_MEAN=1 /D_CRT_SECURE_NO_WARNINGS)
endif()

if(MSVC)
	add_definitions(/bigobj)
else()
	add_definitions(-fwrapv)
endif()
add_definitions(-DUSE_GLOP -DUSE_BOP)

# Verify Dependencies
find_package(Threads REQUIRED)

find_package(Protobuf REQUIRED)
find_package(gflags REQUIRED)
find_package(glog REQUIRED)
find_package(Cbc REQUIRED)

# Main Target
add_library(${PROJECT_NAME} SHARED "")
if(NOT APPLE)
	set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
else()
	# Clang don't support version x.y.z with z > 255
	set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_${PROJECT_NAME}_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPATIBLE_INTERFACE_STRING ${PROJECT_NAME}_MAJOR_VERSION)
if(APPLE)
	set_target_properties(${PROJECT_NAME} PROPERTIES
		INSTALL_RPATH
		"@loader_path")
endif()
target_include_directories(${PROJECT_NAME} INTERFACE
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
	$<INSTALL_INTERFACE:include>
	)
target_link_libraries(${PROJECT_NAME} PUBLIC protobuf::libprotobuf gflags glog Cbc ${CMAKE_THREAD_LIBS_INIT})
target_compile_definitions(${PROJECT_NAME} PUBLIC	USE_BOP USE_GLOP)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# Generate Protobuf cpp sources
set(PROTO_HDRS)
set(PROTO_SRCS)
file(GLOB_RECURSE proto_files RELATIVE ${PROJECT_SOURCE_DIR} "ortools/*.proto")
foreach (PROTO_FILE ${proto_files})
	#message(STATUS "protoc proto: ${PROTO_FILE}")
	get_filename_component(PROTO_DIR ${PROTO_FILE} DIRECTORY)
	get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
	set(PROTO_HDR ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.h)
	set(PROTO_SRC ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.cc)
	#message(STATUS "protoc hdr: ${PROTO_HDR}")
	#message(STATUS "protoc src: ${PROTO_SRC}")
	add_custom_command(
		OUTPUT ${PROTO_SRC} ${PROTO_HDR}
		COMMAND protobuf::protoc
		"--proto_path=${PROJECT_SOURCE_DIR}"
		"--cpp_out=${PROJECT_BINARY_DIR}"
		${PROTO_FILE}
		DEPENDS ${PROTO_FILE} protobuf::protoc
		COMMENT "Running C++ protocol buffer compiler on ${PROTO_FILE}"
		VERBATIM)
	list(APPEND PROTO_HDRS ${PROTO_HDR})
	list(APPEND PROTO_SRCS ${PROTO_SRC})
endforeach()
#add_library(${PROJECT_NAME}_proto STATIC ${PROTO_SRCS} ${PROTO_HDRS})
add_library(${PROJECT_NAME}_proto OBJECT ${PROTO_SRCS} ${PROTO_HDRS})
set_target_properties(${PROJECT_NAME}_proto PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(${PROJECT_NAME}_proto PRIVATE
	${PROJECT_SOURCE_DIR}
	${PROJECT_BINARY_DIR}
	$<TARGET_PROPERTY:protobuf::libprotobuf,INTERFACE_INCLUDE_DIRECTORIES>
	)
#target_link_libraries(${PROJECT_NAME}_proto PRIVATE protobuf::libprotobuf)
add_dependencies(${PROJECT_NAME}_proto protobuf::libprotobuf)
add_library(${PROJECT_NAME}::proto ALIAS ${PROJECT_NAME}_proto)
# Add ortools::proto to libortools
#target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}::proto)
target_sources(${PROJECT_NAME} PRIVATE $<TARGET_OBJECTS:${PROJECT_NAME}::proto>)
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}::proto)

foreach(SUBPROJECT
		algorithms base bop	constraint_solver	data glop	graph	linear_solver	lp_data
		port sat util)
	add_subdirectory(ortools/${SUBPROJECT})
	#target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}::${SUBPROJECT})
	target_sources(${PROJECT_NAME} PRIVATE $<TARGET_OBJECTS:${PROJECT_NAME}::${SUBPROJECT}>)
	add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}::${SUBPROJECT})
endforeach()

# Install rules
include(GNUInstallDirs)

include(GenerateExportHeader)
GENERATE_EXPORT_HEADER(${PROJECT_NAME})
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_export.h
	DESTINATION	${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${PROJECT_NAME}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	)
install(DIRECTORY ortools
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "*.h")
install(DIRECTORY ${PROJECT_BINARY_DIR}/ortools
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "*.pb.h"
	PATTERN CMakeFiles EXCLUDE)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${PROJECT_BINARY_DIR}/ortools/${PROJECT_NAME}ConfigVersion.cmake"
	COMPATIBILITY SameMajorVersion
	)
install(
	FILES
	"${PROJECT_SOURCE_DIR}/ortools/cmake/${PROJECT_NAME}Config.cmake"
	"${PROJECT_BINARY_DIR}/ortools/${PROJECT_NAME}ConfigVersion.cmake"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	COMPONENT Devel)
