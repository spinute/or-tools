set(gflags_URL https://github.com/gflags/gflags)

if (MSVC)
    set(gflags_ADDITIONAL_CMAKE_OPTIONS "-G \"NMake MakeFiles\"")
endif()

set_property(DIRECTORY PROPERTY EP_BASE dependencies)

ExternalProject_Add(gflags_project
	GIT_REPOSITORY ${gflags_URL}
	GIT_TAG "v${gflags_VERSION}"
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>
	-DBUILD_STATIC_LIBS=ON
	-DBUILD_TESTING=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE=ON
	${gflags_ADDITIONAL_CMAKE_OPTIONS}
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

ExternalProject_Get_Property(gflags_project source_dir)
ExternalProject_Get_Property(gflags_project binary_dir)

# Old way
set(gflags_INCLUDE_DIRS ${binary_dir}/include)
set(gflags_LIBRARIES gflags)

# Library
add_library(gflags STATIC IMPORTED)
set_target_properties(gflags PROPERTIES IMPORTED_LOCATION
	${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags.a)
# INTERFACE_INCLUDE_DIRECTORIES does not allow non-existent directories
# cf https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${binary_dir}/include)
set_target_properties(gflags PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
	${binary_dir}/include/)
# Can't Alias imported target.
#add_library(gflags::gflags ALIAS gflags)
add_dependencies(gflags gflags_project)

# Install Rules
include(GNUInstallDirs)
install(FILES
	$<TARGET_PROPERTY:gflags,IMPORTED_LOCATION>
	DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY
	$<TARGET_PROPERTY:gflags,INTERFACE_INCLUDE_DIRECTORIES>/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "*.h"
	PATTERN "CMakeFiles" EXCLUDE
	)
