# setup ROOT includes and libraries
find_package(ROOT)

include_directories(${YAP_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

set(PROGRAMS D4piTest D3piTest DKKpiTest)
foreach(program ${PROGRAMS})
	add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx)
	target_link_libraries(${program} YAP ${ROOT_LIBS})
endforeach(program)
