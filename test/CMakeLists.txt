# setup ROOT includes and libraries
find_package(ROOT)

include_directories(${YAP_SOURCE_DIR}/include ${YAP_SOURCE_DIR}/tools ${ROOT_INCLUDE_DIR})

set(YAP_TEST_SOURCES
  test_ClebschGordan.cxx
  test_FourMomentaCalculation.cxx
  test_HelicityAngles.cxx
#  test_Matrix.cxx
  test_swapDalitzAxes.cxx
  test_Vector.cxx
  test_WignerD.cxx
)

add_executable(yap_test yap_test.cxx ${YAP_TEST_SOURCES})
target_link_libraries(yap_test YAP "${ROOT_LIBS}")


