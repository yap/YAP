cmake_minimum_required(VERSION 2.8.12)
project(YAP)

# Add a custom 'doc' target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

# DataPoint type flag
if (YAP_DATA_POINT_TYPE STREQUAL "float")
	message("Using float as type for DataPoint")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D YAP_DATA_POINT_FLOAT")
else()
	message("Using double as type for DataPoint (to use floats in order to reduce memory usage by DataPoints, cmake with -DYAP_DATA_POINT_TYPE=float)")
endif()

# Default C++ flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
# CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
# elseif(COMPILER_SUPPORTS_CXX0X)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

CHECK_CXX_COMPILER_FLAG("-fdiagnostics-color=always" COMPILER_SUPPORTS_DIAGNOSTIC_COLORS)
if (COMPILER_SUPPORTS_DIAGNOSTIC_COLORS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-long-long -Wno-missing-braces -Werror=overloaded-virtual -pthread")

set(EXT_PROJECTS_DIR  ${YAP_SOURCE_DIR}/external)

# Default build type is 'Release' (cmake default one is '')
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# build with cmake -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g2")

# build with cmake -DCMAKE_BUILD_TYPE=Release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -D ELPP_DISABLE_DEBUG_LOGS")

# Add subdirectories
add_subdirectory(src)
add_subdirectory(test)

# Examples are not part of the `all` target.
# This will provide an `example` target.
add_subdirectory(examples)
