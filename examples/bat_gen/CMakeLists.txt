find_package(ROOT)
find_package(BAT)

include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

set(SOURCES
  bat_fit.cxx
	bat_gen.cxx
  bat_yap_base.cxx
  fit_fitFraction.cxx
	)

#add_custom_target(create-examples-run_gen-output COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/output)

message(STATUS ${CMAKE_CURRENT_BINARY_DIR} )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_custom_target(link-macros
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		../../../examples/bat_gen/macros
		macros
	)

set(PROGRAMS runBatGen runBatFit runBatFitFitFractions runBatFitFitFractionsD4pi)
foreach(program ${PROGRAMS})
	add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx ${SOURCES})
	target_link_libraries(${program} ${ROOT_LIBS} ${BAT_LIBS} YAP)
endforeach(program)

# add_executable(runBatGen EXCLUDE_FROM_ALL runBatGen.cxx ${SOURCES})
# add_dependencies(runBatGen link-macros)
# target_link_libraries(runBatGen ${ROOT_LIBS} ${BAT_LIBS} YAP)
# 
# add_executable(runBatFit EXCLUDE_FROM_ALL runBatGen.cxx ${SOURCES})
# add_dependencies(runBatFit link-macros)
# target_link_libraries(runBatFit ${ROOT_LIBS} ${BAT_LIBS} YAP)
