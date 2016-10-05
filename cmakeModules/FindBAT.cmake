# - Finds BAT instalation
# This module sets up BAT information
# It defines:
# BAT_FOUND          If the BAT is found
# BAT_INCLUDE_DIR    PATH to the include directory
# BAT_INCLUDE_DIRS   PATH to the include directories (not cached)
# BAT_LIBRARIES      BAT libraries
# BAT_LIBRARY_DIR    PATH to the library directory
# BAT_LIBRARY_DIRS   PATH to the library directories (not cached)
#
# Updated by Paolo Di Giglio (p.digiglio91@gmail.com)

find_program(BAT_CONFIG_EXECUTABLE bat-config)

execute_process(
    COMMAND ${BAT_CONFIG_EXECUTABLE} --prefix
    OUTPUT_VARIABLE BATSYS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${BAT_CONFIG_EXECUTABLE} --version
    OUTPUT_VARIABLE BAT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${BAT_CONFIG_EXECUTABLE} --incdir
    OUTPUT_VARIABLE BAT_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${BAT_CONFIG_EXECUTABLE} --libdir
    OUTPUT_VARIABLE BAT_LIBRARY_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${BAT_CONFIG_EXECUTABLE} --libs
    OUTPUT_VARIABLE BAT_LIBRARY
    OUTPUT_STRIP_TRAILING_WHITESPACE)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BAT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BAT DEFAULT_MSG BAT_LIBRARY BAT_INCLUDE_DIR)

mark_as_advanced(BAT_INCLUDE_DIR BAT_LIBRARY)

set(BAT_INCLUDE_DIRS ${BAT_INCLUDE_DIR})
set(BAT_LIBRARY_DIRS ${BAT_INCLUDE_DIR})
set(BAT_LIBRARIES    ${BAT_LIBRARY})
