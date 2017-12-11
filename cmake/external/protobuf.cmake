set(Protobuf_URL https://github.com/google/protobuf)
find_package(ZLIB REQUIRED)

if (MSVC)
    set(Protobuf_ADDITIONAL_CMAKE_OPTIONS "${Protobuf_ADDITIONAL_CMAKE_OPTIONS} -G \"NMake MakeFiles\"")
endif()

set_property(DIRECTORY PROPERTY EP_BASE dependencies)

ExternalProject_Add(Protobuf_project
	GIT_REPOSITORY ${Protobuf_URL}
	GIT_TAG "v${Protobuf_VERSION}"
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>/cmake
	-DBUILD_STATIC_LIBS=ON
	-Dprotobuf_BUILD_TESTS=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE=ON
	${Protobuf_ADDITIONAL_CMAKE_OPTIONS}
	INSTALL_COMMAND ""
	TEST_COMMAND ""
	CMAKE_CACHE_ARGS
	-DCMAKE_BUILD_TYPE:STRING=Release
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
	LOG_DOWNLOAD ON
	LOG_CONFIGURE ON
	LOG_BUILD ON
	)

# Specify include dir
ExternalProject_Get_Property(Protobuf_project source_dir)
ExternalProject_Get_Property(Protobuf_project binary_dir)

# Old way
set(Protobuf_INCLUDE_DIRS ${source_dir}/src)
set(Protobuf_LIBRARIES Protobuf ${ZLIB_LIBRARIES})

# Library
add_library(Protobuf STATIC IMPORTED)
set_target_properties(Protobuf PROPERTIES IMPORTED_LOCATION
	${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}protobuf.a)
# INTERFACE_INCLUDE_DIRECTORIES does not allow non-existent directories
# cf https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${source_dir}/src)
set_target_properties(Protobuf PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
	${source_dir}/src)
set_target_properties(Protobuf PROPERTIES INTERFACE_LINK_LIBRARIES
	${ZLIB_LIBRARIES})
# Can't Alias imported target.
#add_library(Protobuf::Protobuf ALIAS Protobuf)
add_dependencies(Protobuf Protobuf_project)

# Install Rules
include(GNUInstallDirs)
install(FILES
	$<TARGET_PROPERTY:Protobuf,IMPORTED_LOCATION>
	$<TARGET_PROPERTY:Protobuf,INTERFACE_LINK_LIBRARIES>
	DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY
	$<TARGET_PROPERTY:Protobuf,INTERFACE_INCLUDE_DIRECTORIES>/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "google/protobuf/*.h"
	)

# Export protoc executable to consume .proto file during build
add_executable(protobuf::protoc IMPORTED)
set_target_properties(protobuf::protoc PROPERTIES IMPORTED_LOCATION
	"${binary_dir}/protoc")
add_dependencies(protobuf::protoc Protobuf_project)
