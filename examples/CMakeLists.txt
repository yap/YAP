add_subdirectory(programs)
add_subdirectory(bat_gen)

# TODO find a smart way to list all example executables!!
add_custom_target(examples 
	DEPENDS
  	${CMAKE_CURRENT_BINARY_DIR}/bat_gen/runBatGen
  	${CMAKE_CURRENT_BINARY_DIR}/bat_gen/runBatFit
  	${CMAKE_CURRENT_BINARY_DIR}/bat_gen/runBatFitFitFractions
    ${CMAKE_CURRENT_BINARY_DIR}/bat_gen/runBatFitFitFractionsD4pi
	${CMAKE_CURRENT_BINARY_DIR}/programs/D3piTest 
	${CMAKE_CURRENT_BINARY_DIR}/programs/D4piTest 
	${CMAKE_CURRENT_BINARY_DIR}/programs/DKKpiTest 
	)
