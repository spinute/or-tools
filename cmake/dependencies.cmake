enable_language(CXX)

if (MSVC)
    add_definitions(/bigobj /DNDEBUG /DUSE_GLOP /DUSE_BOP)
else()
    add_definitions(-fwrapv -DNDEBUG -DUSE_GLOP -DUSE_BOP)
endif()

if (NOT BUILD_DEPS)
	return()
endif()

# Find or build Dependencies
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Dependencies.txt _Dependency_file)
foreach(DEPENDENCY Protobuf gflags glog Cbc CoinUtils Osi Clp Cgl)
    string(REGEX REPLACE ".*${DEPENDENCY} = ([0-9.]+).*" "\\1" ${DEPENDENCY}_VERSION ${_Dependency_file})
endforeach()

if (MSVC)
	include(external/zlib)
	include(external/swig)
endif()

# Protobuf
if (BUILD_DEPS)
	set(Protobuf_FOUND False)
else()
	find_package(Protobuf ${Protobuf_VERSION})
endif()
if (NOT Protobuf_FOUND)
	message(STATUS "Did not find system protobuf or forced build. Building as an external project")
	include(external/protobuf)
endif()
#include_directories(${Protobuf_INCLUDE_DIRS})

# gflags
if (BUILD_DEPS)
	set(gflags_FOUND False)
else()
	find_package(gflags ${gflags_VERSION})
endif()
if (NOT gflags_FOUND)
	message(STATUS "Did not find system gflags or forced build. Building as an external project")
	include(external/gflags)
endif()
#include_directories(${gflags_INCLUDE_DIRS})

# glog
if (BUILD_DEPS)
	set(glog_FOUND False)
else()
	find_package(glog ${glog_VERSION})
endif()
if (NOT glog_FOUND)
	message(STATUS "Did not find system glog or forced build. Building as an external project.")
	include(external/glog)
endif()
#include_directories(${glog_INCLUDE_DIRS})

# Cbc
if (BUILD_DEPS)
	set(Cbc_FOUND False)
else()
	find_package(Cbc ${Cbc_VERSION})
endif()
if (NOT Cbc_FOUND)
	#if (NOT MSVC)
	message(STATUS "Did not find system Cbc or forced build. Building as an external project.")
	include(external/cbc)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CLP -DUSE_CBC")
	#include_directories(${Cbc_INCLUDE_DIRS})
	#endif()
	#else()
endif()
#include_directories(${Cbc_INCLUDE_DIRS})
#if (MSVC)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DUSE_CLP /DUSE_CBC")
#else()
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CLP -DUSE_CBC")
#endif()
#add_definitions("-DUSE_CLP -DUSE_CBC")
