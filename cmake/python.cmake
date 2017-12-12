if (NOT BUILD_PYTHON)
	return()
endif()

if(NOT TARGET ortools::ortools)
	message(FATAL_ERROR "Python: missing ortools TARGET")
endif()

# Will need swig
find_package(SWIG REQUIRED)
include(UseSWIG)

# Generate Protobuf py sources
set(PROTO_PYS)
file(GLOB_RECURSE proto_py_files RELATIVE ${PROJECT_SOURCE_DIR}
	"ortools/constraint_solver/*.proto"
	"ortools/linear_solver/*.proto")
list(REMOVE_ITEM proto_py_files "ortools/constraint_solver/demon_profiler.proto")
foreach(PROTO_FILE ${proto_py_files})
	message(STATUS "protoc: ${PROTO_FILE}")
	get_filename_component(PROTO_DIR ${PROTO_FILE} DIRECTORY)
	get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
	set(PROTO_PY ${PROJECT_BINARY_DIR}/${PROTO_DIR}/${PROTO_NAME}_pb2.py)
	message(STATUS "protoc py: ${PROTO_PY}")
	add_custom_command(
		OUTPUT ${PROTO_PY}
		COMMAND protobuf::protoc
		"--proto_path=${PROJECT_SOURCE_DIR}"
		"--python_out=${PROJECT_BINARY_DIR}"
		${PROTO_FILE}
		DEPENDS ${PROTO_FILE} protobuf::protoc
		COMMENT "Running C++ protocol buffer compiler on ${PROTO_FILE}"
		VERBATIM)
	list(APPEND PROTO_PYS ${PROTO_PY})
endforeach()
add_custom_target(Py${PROJECT_NAME}_proto DEPENDS ${PROTO_PYS} ortools::ortools)

# Build Python wrappers
# Specify supported python version
set(Python_ADDITIONAL_VERSIONS "3.6;3.5;2.7"
	CACHE STRING "available python version")
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

if(${PYTHON_VERSION_STRING} VERSION_GREATER 3)
	set(CMAKE_SWIG_FLAGS "-py3;-DPY3")
endif()
list(APPEND CMAKE_SWIG_FLAGS "-I${PROJECT_SOURCE_DIR}")

foreach(SUBPROJECT constraint_solver linear_solver sat graph algorithms data)
	add_subdirectory(ortools/${SUBPROJECT}/python)
endforeach()

# Copy to binary dir
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/)
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/constraint_solver/)
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/linear_solver/)
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/sat/)
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/graph/)
file(COPY
	${PROJECT_SOURCE_DIR}/ortools/__init__.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/algorithms/)

file(COPY
	${PROJECT_SOURCE_DIR}/ortools
	DESTINATION
	${PROJECT_BINARY_DIR}
	FILES_MATCHING
	PATTERN
	"*.i")

file(COPY
	${PROJECT_SOURCE_DIR}/README
	DESTINATION
	${PROJECT_BINARY_DIR})

file(COPY
	${PROJECT_SOURCE_DIR}/ortools/linear_solver/linear_solver_natural_api.py
	DESTINATION
	${PROJECT_BINARY_DIR}/ortools/linear_solver/)

file(COPY
	${PROJECT_SOURCE_DIR}/python/MANIFEST.in
	DESTINATION
	${PROJECT_BINARY_DIR}/)
set(README_FILE README)

# Main Target
configure_file(${PROJECT_SOURCE_DIR}/python/setup.py.in ${PROJECT_BINARY_DIR}/setup.py)
set(PY_OUTPUT ${PROJECT_BINARY_DIR}/timestamp)
add_custom_command(
	OUTPUT ${PY_OUTPUT}
	COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/setup.py sdist
	COMMAND ${CMAKE_COMMAND} -E touch ${PY_OUTPUT}
	DEPENDS Py${PROJECT_NAME}_proto ${PROJECT_SOURCE_DIR}/ortools/__init__.py)
add_custom_target(Py${PROJECT_NAME} ALL
	DEPENDS ${PY_OUTPUT})

# Install rules
include(GNUInstallDirs)

