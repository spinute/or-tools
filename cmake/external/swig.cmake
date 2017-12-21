set(SWIG_URL http://prdownloads.sourceforge.net/swig/swigwin-3.0.12.zip)

set_property(DIRECTORY PROPERTY EP_BASE dependencies)

ExternalProject_Add(SWIG_project
	URL ${SWIG_URL}
	BUILD_IN_SOURCE 1
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	TEST_COMMAND ""
	LOG_DOWNLOAD ON
	LOG_CONFIGURE ON
	LOG_BUILD ON
	)

# Specify include dir
ExternalProject_Get_Property(SWIG_project source_dir)
ExternalProject_Get_Property(SWIG_project binary_dir)

# Executable
add_executable(SWIG_EXECUTABLE STATIC IMPORTED)
set_target_properties(SWIG_EXECUTABLE PROPERTIES IMPORTED_LOCATION
	${binary_dir}/swigwin-3.0.12/swig.exe)
add_dependencies(SWIG_EXECUTABLE SWIG_project)

