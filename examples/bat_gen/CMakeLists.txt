find_package(ROOT)
find_package(BAT )

include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

set(SOURCES
	runBatGen.cxx
	bat_gen.cxx
	)

add_executable(runBatGen EXCLUDE_FROM_ALL ${SOURCES})
target_link_libraries(runBatGen ${ROOT_LIBS} ${BAT_LIBS} YAP)
