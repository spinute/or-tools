if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE sat.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapsat")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapsat TYPE SHARED LANGUAGE python SOURCES sat.i)
else()
	swig_add_module(pywrapsat python sat.i)
endif()
swig_link_libraries(pywrapsat ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapsat_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})

add_dependencies(Py${PROJECT_NAME}_proto ${SWIG_MODULE_pywrapsat_REAL_NAME})

#[[
include(GNUInstallDirs)
install(TARGETS ${SWIG_MODULE_pywrapsat_REAL_NAME}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
]]
