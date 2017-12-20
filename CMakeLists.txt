cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
project(ortools-meta NONE)

include(CMakeDependentOption)
option(BUILD_DEPS "Force re-build of all dependencies" ON)
set(DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/download" CACHE PATH "Location where external projects will be downloaded.")
mark_as_advanced(DOWNLOAD_LOCATION)

option(BUILD_CXX "Build C++ library" ON)
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON "Build Python Library" ON "BUILD_CXX"	OFF)
CMAKE_DEPENDENT_OPTION(BUILD_JAVA "Build Java Library" OFF "BUILD_CXX" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_CSHARP "Build CSharp Library" OFF "BUILD_CXX" OFF)

message(STATUS "Build all dependencies: ${BUILD_DEPS}")
message(STATUS "Build CXX library: ${BUILD_CXX}")
message(STATUS "Build Python library: ${BUILD_PYTHON}")
message(STATUS "Build Java library: ${BUILD_JAVA}")
message(STATUS "Build CSharp library: ${BUILD_CSHARP}")

include(CTest)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Required dependencies
include(dependencies)

cmake_policy(SET CMP0068 NEW)
#add_subdirectory(ortools)
include(cpp)
include(python)
include(java)
include(csharp)

add_subdirectory(examples/cpp)
add_subdirectory(examples/data)
add_subdirectory(examples/python)
add_subdirectory(examples/notebook)
