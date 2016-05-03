cmake_minimum_required(VERSION 2.8.12)
project(YAP)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakeModules)

# Default C++ flags
# TODO -std=c++11 is supported only since gcc 4.7: check compiler version!
set(CMAKE_CXX_FLAGS "-std=c++11 -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-long-long -Wno-missing-braces -Werror=overloaded-virtual -pthread" )

set(EXT_PROJECTS_DIR  ${YAP_SOURCE_DIR}/external)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Building ${CMAKE_BUILD_TYPE} project")

# build with cmake -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g2" )

# build with cmake -DCMAKE_BUILD_TYPE=Release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -D ELPP_DISABLE_DEBUG_LOGS" )

# Add subdirectories
add_subdirectory(src)
add_subdirectory(examples)
#add_subdirectory(test)
