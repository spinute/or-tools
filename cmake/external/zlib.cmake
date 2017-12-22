# Download and unpack zlib at configure time
configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeLists.txt.zlib zlib-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/zlib-download )
if(result)
  message(FATAL_ERROR "CMake step for zlib failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/zlib-download )
if(result)
  message(FATAL_ERROR "Build step for zlib failed: ${result}")
endif()

# Old way
set(ZLIB_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/zlib/)
set(ZLIB_LIBRARIES ZLIB)

# Library
add_library(ZLIB STATIC IMPORTED)
set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION
	${CMAKE_BINARY_DIR}/zlib/zlib.lib)
set_target_properties(ZLIB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
	${CMAKE_BINARY_DIR}/zlib)
# Can't Alias imported target.
#add_library(ZLIB::ZLIB ALIAS ZLIB)

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

