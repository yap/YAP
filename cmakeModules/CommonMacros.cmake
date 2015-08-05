function(make_executable EXE_NAME SOURCES)
	if(DEBUG_OUTPUT)
		message(STATUS "make_executable was called with the following arguments:
        EXE_NAME = '${EXE_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
	endif()
	add_executable(${EXE_NAME} ${SOURCES})
	# proccess link libraries in additional arguments
	foreach(_LIB ${ARGN})
		target_link_libraries(${EXE_NAME} ${_LIB})
	endforeach()
	unset(_LIB)
endfunction(make_executable)

function(make_shared_library LIB_NAME SOURCES)
	if(DEBUG_OUTPUT)
		message(STATUS "make_shared_library was called with the following arguments:
        LIB_NAME = '${LIB_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
	endif()
	add_library(${LIB_NAME} SHARED ${SOURCES})
	# proccess link libraries in additional arguments
	foreach(_LIB ${ARGN})
		target_link_libraries(${LIB_NAME} ${_LIB})
	endforeach()
	unset(_LIB)
endfunction(make_shared_library)

function(make_static_library LIB_NAME SOURCES)
	if(DEBUG_OUTPUT)
		message(STATUS "make_static_library was called with the following arguments:
        LIB_NAME = '${LIB_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
	endif()
	add_library(${LIB_NAME} STATIC ${SOURCES})
	# proccess link libraries in additional arguments
	foreach(_LIB ${ARGN})
		target_link_libraries(${LIB_NAME} ${_LIB})
	endforeach()
	unset(_LIB)
endfunction(make_static_library)