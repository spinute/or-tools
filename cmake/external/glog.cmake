set(glog_URL https://github.com/google/glog)

ExternalProject_Get_Property(gflags_project source_dir)
set(glog_ADDITIONAL_CMAKE_OPTIONS "-DCMAKE_PREFIX_PATH=${source_dir}")

if (MSVC)
    set(glog_ADDITIONAL_CMAKE_OPTIONS "-G \"NMake MakeFiles\"")
endif()

set_property(DIRECTORY PROPERTY EP_BASE dependencies)

ExternalProject_Add(glog_project
	GIT_REPOSITORY ${glog_URL}
	GIT_TAG "v${glog_VERSION}"
	UPDATE_COMMAND ""
	PATCH_COMMAND git am -3 ${CMAKE_SOURCE_DIR}/patches/glog_includedir_fix.patch
	CONFIGURE_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR>
	-DWITH_GFLAGS=ON
	-DBUILD_TESTING=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE=ON
	${glog_ADDITIONAL_CMAKE_OPTIONS}
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

add_dependencies(glog_project gflags_project)

ExternalProject_Get_Property(glog_project source_dir)
ExternalProject_Get_Property(glog_project binary_dir)

# Old way
set(glog_INCLUDE_DIRS ${source_dir}/src ${binary_dir})
set(glog_LIBRARIES glog)

# Library
add_library(glog STATIC IMPORTED)
set_target_properties(glog PROPERTIES IMPORTED_LOCATION
	${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}glog.a)
# INTERFACE_INCLUDE_DIRECTORIES does not allow non-existent directories
# cf https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${source_dir}/src ${binary_dir})
set_target_properties(glog PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
	"${source_dir}/src/;${binary_dir}/")
# Can't Alias imported target.
#add_library(glog::glog ALIAS glog)
add_dependencies(glog glog_project)

# Install Rules
include(GNUInstallDirs)
install(FILES
	$<TARGET_PROPERTY:glog,IMPORTED_LOCATION>
	DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY
	$<TARGET_PROPERTY:glog,INTERFACE_INCLUDE_DIRECTORIES>
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "*.h"
	PATTERN "CMakeFiles" EXCLUDE
	)
