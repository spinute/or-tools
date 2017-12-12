cmake_minimum_required(VERSION 3.5)
project(ortools-meta NONE)

include(CTest)
set(CMAKE_VERBOSE_MAKEFILE ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(CMakeDependentOption)

option(BUILD_DEPS "Force re-build of all dependencies" ON)
set(DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/download" CACHE PATH "Location where external projects will be downloaded.")
mark_as_advanced(DOWNLOAD_LOCATION)

option(BUILD_CXX "Build C++ library" ON)
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON "Build Python Library" OFF "NOT BUILD_CXX" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_JAVA "Build Java Library" OFF "NOT BUILD_CXX" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_CSHARP "Build CSharp Library" OFF "NOT BUILD_CXX" OFF)

message(STATUS "Build all dependencies: ${BUILD_DEPS}")
message(STATUS "Build CXX library: ${BUILD_CXX}")
message(STATUS "Build Python library: ${BUILD_PYTHON}")
message(STATUS "Build Java library: ${BUILD_JAVA}")
message(STATUS "Build CSharp library: ${BUILD_CSHARP}")

# For Debug only ?
set(CMAKE_INSTALL_PREFIX install
	CACHE STRING "install prefix for debug purpose")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	CACHE STRING "Build directory for .a")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	CACHE STRING "Build directory for .so")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
	CACHE STRING "Build directory for binary")

# Required dependencies
include(dependencies)

#add_subdirectory(ortools)
include(cpp)
include(python)
include(java)
include(csharp)

add_subdirectory(examples/cpp)
add_subdirectory(examples/data)
add_subdirectory(examples/python)
