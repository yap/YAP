cmake_minimum_required(VERSION 2.8.8)
project(YAP)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakeModules)
set(CMAKE_CXX_FLAGS "-std=c++11" )

include(CommonMacros)

# setup ROOT includes and libraries
find_package(ROOT 5.34)

#BAT
# find_program(BAT_CONFIG_EXECUTABLE bat-config)
# if (BAT_CONFIG_EXECUTABLE)
# 	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --incdir
# 		OUTPUT_VARIABLE BAT_INCLUDE_DIR
# 		OUTPUT_STRIP_TRAILING_WHITESPACE)
# 
# 	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --libdir
# 		OUTPUT_VARIABLE BAT_LIBRARY_DIR
# 		OUTPUT_STRIP_TRAILING_WHITESPACE)
# 
# 	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -I${BAT_INCLUDE_DIR})
# 	set(LIBS ${LIBS} -L${BAT_LIBRARY_DIR} -lBAT)
# endif()

if (CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
    list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
    list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

    message(STATUS ${GCC_MAJOR})
    message(STATUS ${GCC_MINOR})
endif()
#set( CMAKE_CXX_COMPILER /nfs/mnemosyne/sys/slc6/contrib/gcc/4.9.0/bin/c++)
#set( CMAKE_CXX_COMPILER clang)
set( CMAKE_CXX_COMPILER g++)
message(${CMAKE_CXX_COMPILER})

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-long-long -Werror=overloaded-virtual" )

add_subdirectory(src)
add_subdirectory(programs)

