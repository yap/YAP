set(BAT_FOUND        FALSE)
set(BAT_ERROR_REASON "")
set(BAT_DEFINITIONS  "")
set(BAT_LIBS)


find_program(BAT_CONFIG_EXECUTABLE bat-config)
if(NOT BAT_CONFIG_EXECUTABLE)
	set(BAT_ERROR_REASON "${BAT_ERROR_REASON} $Cannot find bat-config executable in path. Make sure BAT is setup correctly.")
else()

	set(BAT_FOUND TRUE)

	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --prefix
		OUTPUT_VARIABLE BATSYS
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --version
		OUTPUT_VARIABLE BAT_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --incdir
		OUTPUT_VARIABLE BAT_INCLUDE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT EXISTS "${BAT_INCLUDE_DIR}")
		set(BAT_FOUND FALSE)
		set(BAT_ERROR_REASON "${BAT_ERROR_REASON} BAT include directory '${BAT_INCLUDE_DIR}' does not exist.")
	endif()

	execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --libdir
		OUTPUT_VARIABLE BAT_LIBRARY_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT EXISTS "${BAT_LIBRARY_DIR}")
		set(BAT_FOUND FALSE)
		set(BAT_ERROR_REASON "${BAT_ERROR_REASON} BAT library directory '${BAT_LIBRARY_DIR}' does not exist.")
		set(BAT_LIBS "${BAT_LIBRARY_DIR} -lBAT")
	endif()

endif()

# make variables changeable
mark_as_advanced(
	BAT_INCLUDE_DIR
	BAT_LIBRARY_DIR
	BAT_LIBRARIES
	)


# report result
if(BAT_FOUND)
	message(STATUS "Found BAT version ${BAT_VERSION} r${BAT_SVN_REVISION} in '${BATSYS}'.")
	message(STATUS "Using BAT include directory '${BAT_INCLUDE_DIR}'.")
	message(STATUS "Using BAT library directory '${BAT_LIBRARY_DIR}'.")
	message(STATUS "Using BAT libraries ${BAT_LIBRARIES}.")
else()
		message(FATAL_ERROR "Unable to find requested BAT installation:${BAT_ERROR_REASON}")
endif()