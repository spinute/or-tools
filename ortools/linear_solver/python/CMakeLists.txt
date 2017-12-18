if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE linear_solver.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywraplp")

# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/..)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../../:$ORIGIN/")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywraplp TYPE SHARED LANGUAGE python SOURCES linear_solver.i)
else()
	swig_add_module(pywraplp python linear_solver.i)
endif()
swig_link_libraries(pywraplp ortools::ortools)
target_include_directories(${SWIG_MODULE_pywraplp_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
