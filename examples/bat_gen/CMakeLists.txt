#find_package(ROOT) # Done in the ${YAPDIR}/CMakeLists.txt
#find_package(BAT)  #  ""

include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

set(SOURCES
	runBatGen.cxx
	bat_gen.cxx
	)

message(STATUS "ROOT_INCLUDE_DIR ${ROOT_INCLUDE_DIR}")
add_executable(runBatGen ${SOURCES})
message(STATUS "${ROOT_LIBRARIES} ${BAT_LIBS} YAP")
target_link_libraries(runBatGen ${ROOT_LIBS} ${BAT_LIBS} YAP)
