if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE linear_solver.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywraplp")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywraplp TYPE SHARED LANGUAGE python SOURCES linear_solver.i)
else()
	swig_add_module(pywraplp python linear_solver.i)
endif()
swig_link_libraries(pywraplp ortools::ortools)
target_include_directories(${SWIG_MODULE_pywraplp_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})

add_dependencies(Py${PROJECT_NAME}_proto ${SWIG_MODULE_pywraplp_REAL_NAME})

#[[
include(GNUInstallDirs)
install(TARGETS ${SWIG_MODULE_pywraplp_REAL_NAME}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
]]
