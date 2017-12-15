if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE knapsack_solver.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywrapknapsack_solver")

# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/..)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywrapknapsack_solver TYPE SHARED LANGUAGE python SOURCES knapsack_solver.i)
else()
	swig_add_module(pywrapknapsack_solver python knapsack_solver.i)
endif()
swig_link_libraries(pywrapknapsack_solver ortools::ortools)
target_include_directories(${SWIG_MODULE_pywrapknapsack_solver_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
