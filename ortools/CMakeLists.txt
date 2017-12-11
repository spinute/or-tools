if (NOT BUILD_CXX)
	return()
endif()

cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(ORToolsHelpers)
set_version(ORTOOLS_VERSION)
project(ortools LANGUAGES CXX VERSION ${ORTOOLS_VERSION})

# SWIG Python Java tools
if (BUILD_PYTHON OR BUILD_JAVA)
	find_package(SWIG REQUIRED)
	include(UseSWIG)
endif()

if (BUILD_PYTHON)
	# Specify python version
	set(Python_ADDITIONAL_VERSIONS "3.6;3.5;2.7"
		CACHE STRING "available python version")
	find_package(PythonInterp REQUIRED)
	find_package(PythonLibs REQUIRED)
endif()

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

# find dependencies
## pthread
find_package(Threads REQUIRED)

check_target(Protobuf)
check_target(gflags)
check_target(glog)
check_target(Cbc)
#check_target(absl)

# Main Target
add_library(${PROJECT_NAME} SHARED "")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_${PROJECT_NAME}_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPATIBLE_INTERFACE_STRING ${PROJECT_NAME}_MAJOR_VERSION)
target_include_directories(${PROJECT_NAME} INTERFACE
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
	)
target_link_libraries(${PROJECT_NAME} PUBLIC Protobuf gflags glog Cbc ${CMAKE_THREAD_LIBS_INIT})

# Project use ortools as root include dir aka #include "ortools/x/y.h" not #include "x/y.h"
get_filename_component(PARENT_SOURCE_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
get_filename_component(PARENT_BINARY_DIR ${PROJECT_BINARY_DIR} DIRECTORY)

# Generate/Build Protobuf sources
set(PROTO_HDRS)
set(PROTO_SRCS)
file(GLOB_RECURSE PROTO_FILES RELATIVE ${PROJECT_SOURCE_DIR} "*.proto")
foreach (PROTO_FILE ${proto_files})
	message(STATUS "protoc proto: ${PROTO_FILE}")
	file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/)
	get_filename_component(PROTO_DIR ${PROTO_FILE} DIRECTORY)
	get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)

	set(PROTO_HDR ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.h)
	set(PROTO_SRC ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}.pb.cc)
	message(STATUS "protoc hdr: ${PROTO_HDR}")
	message(STATUS "protoc src: ${PROTO_SRC}")

	add_custom_command(
		OUTPUT ${PROTO_SRC} ${PROTO_HDR}
		COMMAND protobuf::protoc
		"--cpp_out=${PROJECT_BINARY_DIR}"
		"--proto_path=${PARENT_SOURCE_DIR}"
		${PROTO_FILE}
		DEPENDS ${PROTO_FILE} protobuf::protoc
		COMMENT "Running C++ protocol buffer compiler on ${PROTO}"
		VERBATIM)
	list(APPEND PROTO_HDRS ${PROTO_HDR})
	list(APPEND PROTO_SRCS ${PROTO_SRC})
endforeach()
add_library(${PROJECT_NAME}_proto STATIC ${PROTO_SRCS} ${PROTO_HDRS})
set_target_properties(${PROJECT_NAME}_proto PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CMAKE_CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME}_proto PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
target_include_directories(${PROJECT_NAME}_proto PRIVATE
	${PROJECT_SOURCE_DIR}
	${PROJECT_BINARY_DIR}
	#$<TARGET_PROPERTY:Protobuf,INTERFACE_INCLUDE_DIRECTORIES>
	)
target_link_libraries(${PROJECT_NAME}_proto PRIVATE Protobuf)
add_library(${PROJECT_NAME}::proto ALIAS ${PROJECT_NAME}_proto)

# Add ortools::proto to libortools
target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}::proto)
#target_sources(${PROJECT_NAME} PRIVATE	$<TARGET_OBJECTS:${PROJECT_NAME}::proto>)

foreach(SUBPROJECT
		algorithms base bop	constraint_solver	data glop	graph	linear_solver	lp_data
		port sat util)
	  add_subdirectory(${SUBPROJECT})
		target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}::${SUBPROJECT})
		#target_sources(${PROJECT_NAME} PRIVATE
		#	$<TARGET_OBJECTS:${PROJECT_NAME}::${SUBPROJECT}>)
endforeach()

# Install rules
include(GNUInstallDirs)
include(GenerateExportHeader)
GENERATE_EXPORT_HEADER(${PROJECT_NAME})
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_export.h
	DESTINATION	${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${PROJECT_NAME}
	EXPORT ${PROJECT_NAME}Targets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

install(DIRECTORY ortools
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Devel
    FILES_MATCHING
		PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ortools
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Devel
    FILES_MATCHING
		PATTERN "*.pb.h"
    PATTERN CMakeFiles EXCLUDE)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
	COMPATIBILITY SameMajorVersion
	)
install(
	EXPORT ${PROJECT_NAME}Targets
	NAMESPACE ${PROJECT_NAME}::
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	)
install(
	FILES
	"${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
	"${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	COMPONENT Devel)
