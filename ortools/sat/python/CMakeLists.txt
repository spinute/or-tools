if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE sat.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapsat")

# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/..)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapsat TYPE SHARED LANGUAGE python SOURCES sat.i)
else()
	swig_add_module(pywrapsat python sat.i)
endif()
swig_link_libraries(pywrapsat ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapsat_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
