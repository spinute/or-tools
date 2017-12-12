if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE graph.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapgraph")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapgraph TYPE SHARED LANGUAGE python SOURCES graph.i)
else()
	swig_add_module(pywrapgraph python graph.i)
endif()
swig_link_libraries(pywrapgraph ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapgraph_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})

add_dependencies(Py${PROJECT_NAME}_proto ${SWIG_MODULE_pywrapgraph_REAL_NAME})

#[[
include(GNUInstallDirs)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pywrapgraph.py
#	DESTINATION	${CMAKE_INSTALL_LIBDIR}/)
install(TARGETS ${SWIG_MODULE_pywrapgraph_REAL_NAME}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
]]
