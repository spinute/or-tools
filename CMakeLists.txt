# This file is just an orchestration
cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
		"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)
# Use find_package everywhere, no-op if it's a target
macro(find_package)
	if(NOT TARGET "${ARGV0}")
		_find_package(${ARGV})
	else()
		get_target_property(TGT_VER ${ARGV0} VERSION)
		message(STATUS "Found ${ARGV0}: CMake Target (found version \"${TGT_VER}\")")
		set(${ARGV0}_FOUND TRUE)
	endif()
endmacro()
# Apple: Don't modify install_name when touching RPATH.
if(POLICY CMP0068)
	cmake_policy(SET CMP0068 NEW)
endif()

project(ortools-meta NONE)

include(CTest)
set(CMAKE_VERBOSE_MAKEFILE ON)
if(UNIX)
  # Needed to create python packagefrom the build directory
	set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
	option(BUILD_SHARED_LIBS "Build shared libraries(.so)." ON)
else()
	# Windows only support static build.
 set(BUILD_SHARED_LIBS OFF)
endif()

include(CMakeDependentOption)
option(BUILD_DEPS "Force re-build of all dependencies" ON)
option(BUILD_CXX "Build C++ library" ON)
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON "Build Python Library" ON "BUILD_CXX"	OFF)
CMAKE_DEPENDENT_OPTION(BUILD_JAVA "Build Java Library" OFF "BUILD_CXX" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_CSHARP "Build CSharp Library" OFF "BUILD_CXX" OFF)

message(STATUS "Build all dependencies: ${BUILD_DEPS}")
message(STATUS "Build CXX library: ${BUILD_CXX}")
message(STATUS "Build Python Binding: ${BUILD_PYTHON}")
message(STATUS "Build Java Binding: ${BUILD_JAVA}")
message(STATUS "Build CSharp Binding: ${BUILD_CSHARP}")

# Build Dependencies
add_subdirectory(cmake/external)

set(CMAKE_SKIP_INSTALL_RULES TRUE)

include(cpp)
include(python)
include(java)
include(csharp)

if(ENABLE_TESTING)
	add_subdirectory(examples/cpp)
	add_subdirectory(examples/data)
	add_subdirectory(examples/python)
	add_subdirectory(examples/notebook)
endif()
