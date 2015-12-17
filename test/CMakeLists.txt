include_directories(${YAP_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

set(YAP_TEST_SOURCES
  test_Vector.cxx
	test_WignerD.cxx
)

add_executable(yap_test yap_test.cxx ${YAP_TEST_SOURCES})
target_link_libraries(yap_test YAP)


