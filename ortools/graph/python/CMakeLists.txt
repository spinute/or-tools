if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE graph.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapgraph")

# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/..)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../../:$ORIGIN/")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapgraph TYPE SHARED LANGUAGE python SOURCES graph.i)
else()
	swig_add_module(pywrapgraph python graph.i)
endif()
swig_link_libraries(pywrapgraph ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapgraph_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
