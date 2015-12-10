include_directories(${YAP_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

#set(YAP_TEST_SOURCES
#	catch_WignerD.cxx
#)

#make_static_library(YAP_TEST "${YAP_TEST_SOURCES}")

make_executable(yap_test yap_test YAP "${ROOT_LIBS}")


