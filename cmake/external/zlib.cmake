set(ZLIB_URL https://zlib.net/zlib1211.zip)

set_property(DIRECTORY PROPERTY EP_BASE dependencies)

ExternalProject_Add(ZLIB_project
	URL ${ZLIB_URL}
	BUILD_IN_SOURCE 1
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND "nmake -f zlib-1.2.11\\win32\\Makefile.msc zlib.lib"
	INSTALL_COMMAND ""
	TEST_COMMAND ""
	LOG_DOWNLOAD ON
	LOG_CONFIGURE ON
	LOG_BUILD ON
	)

# Specify include dir
ExternalProject_Get_Property(ZLIB_project source_dir)
ExternalProject_Get_Property(ZLIB_project binary_dir)

# Old way
set(ZLIB_INCLUDE_DIRS ${source_dir}/zlib-1.2.11)
set(ZLIB_LIBRARIES ZLIB)

# Library
add_library(ZLIB STATIC IMPORTED)
set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION
	${binary_dir}/zlib.lib)
# INTERFACE_INCLUDE_DIRECTORIES does not allow non-existent directories
# cf https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${source_dir}/zlib-1.2.11)
set_target_properties(ZLIB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
	${source_dir}/zlib-1.2.11)
# Can't Alias imported target.
#add_library(ZLIB::ZLIB ALIAS ZLIB)
add_dependencies(ZLIB ZLIB_project)

# Install Rules
include(GNUInstallDirs)
install(FILES
	$<TARGET_PROPERTY:ZLIB,IMPORTED_LOCATION>
	DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY
	$<TARGET_PROPERTY:ZLIB,INTERFACE_INCLUDE_DIRECTORIES>
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
	FILES_MATCHING
	PATTERN "*/zlib.h"
	PATTERN "*/zconf.h"
	)

