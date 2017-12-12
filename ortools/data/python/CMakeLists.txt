if(NOT BUILD_PYTHON)
	return()
endif()

set_property(SOURCE rcpsp.i PROPERTY CPLUSPLUS ON)
list(APPEND CMAKE_SWIG_FLAGS "-module;pywraprcpsp")

if(${CMAKE_VERSION} VERSION_GREATER 3.8.2)
	swig_add_library(pywraprcpsp TYPE SHARED LANGUAGE python SOURCES rcpsp.i)
else()
	swig_add_module(pywraprcpsp python rcpsp.i)
endif()
swig_link_libraries(pywraprcpsp ortools::ortools)
target_include_directories(${SWIG_MODULE_pywraprcpsp_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})

add_dependencies(Py${PROJECT_NAME}_proto ${SWIG_MODULE_pywraprcpsp_REAL_NAME})

#message(STATUS "rcpsp: ${swig_outdir}/pywraprcpsp.py")
#$<TARGET_FILE:${SWIG_MODULE_pywraprcpsp_REAL_NAME}>

#[[
include(GNUInstallDirs)
install(TARGETS ${SWIG_MODULE_pywraprcpsp_REAL_NAME}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
]]
