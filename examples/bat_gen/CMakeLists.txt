find_package(ROOT)
find_package(BAT )

include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

set(SOURCES
	runBatGen.cxx
	bat_gen.cxx
	)

#add_custom_target(create-examples-run_gen-output COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/output)

message(STATUS ${CMAKE_CURRENT_BINARY_DIR} )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_custom_target(link-macros
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		../../../examples/bat_gen/macros
		macros
	)

add_executable(runBatGen EXCLUDE_FROM_ALL ${SOURCES})
add_dependencies(runBatGen link-macros)
target_link_libraries(runBatGen ${ROOT_LIBS} ${BAT_LIBS} YAP)
