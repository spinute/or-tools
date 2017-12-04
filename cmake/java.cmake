if (NOT BUILD_JAVA)
	return()
endif()

find_package(SWIG)
find_package(JAVA)
find_package(JNI)

