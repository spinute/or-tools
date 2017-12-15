if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE routing.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapcp")

# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/..)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapcp TYPE SHARED LANGUAGE python SOURCES routing.i)
else()
	swig_add_module(pywrapcp python routing.i)
endif()
swig_link_libraries(pywrapcp ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapcp_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
