include_directories(${YAP_SOURCE_DIR}/include)

set(PROGRAMS 
	D4piTest 
	D3piTest 
	DKKpiTest)
foreach(program ${PROGRAMS})
	add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx)
	target_link_libraries(${program} YAP)
endforeach(program)
