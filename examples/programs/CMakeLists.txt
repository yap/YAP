include_directories(${YAP_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

set(PROGRAMS 
	D4piTest 
	D3piTest 
	DKKpiTest 
	HelicityZemachRelativeScaling 
	HelicityZemachRelativeScaling2
	HelicityZemachRelativeScaling3)
foreach(program ${PROGRAMS})
	add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx)
	target_link_libraries(${program} ${ROOT_LIBS} YAP)
endforeach(program)
