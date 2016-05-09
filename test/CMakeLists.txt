#find_package(ROOT) # Checked for in ${YAPDIR}/CMakeLists.txt
find_package(Git)

include(ExternalProject)
set(CATCH_SOURCE_DIR ${YAP_SOURCE_DIR}/external/catch)
ExternalProject_Add(
	catch
	DOWNLOAD_DIR      ${CATCH_SOURCE_DIR}
	GIT_REPOSITORY    https://github.com/philsquared/Catch
	UPDATE_COMMAND    ${GIT_EXECUTABLE} pull origin master
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
	LOG_DOWNLOAD      ON
)
set_target_properties(catch PROPERTIES EXCLUDE_FROM_ALL TRUE)
ExternalProject_Get_Property(catch source_dir)

include_directories(${YAP_SOURCE_DIR}/include ${source_dir}/include ${YAP_SOURCE_DIR}/tools ${ROOT_INCLUDE_DIR})
set(YAP_TEST_SOURCES
  test_ClebschGordan.cxx
  test_FourMomentaCalculation.cxx
  test_helicityFrame.cxx
  test_HelicityAngles.cxx
  test_HelicityAngles_boostRotate.cxx
  test_Matrix.cxx
  test_swapDalitzAxes.cxx
  test_swapFinalStates.cxx
  test_Vector.cxx
  test_WignerD.cxx
  yap_test.cxx
)

add_executable(yap_test EXCLUDE_FROM_ALL ${YAP_TEST_SOURCES})
add_dependencies(yap_test catch)
target_link_libraries(yap_test YAP ${ROOT_LIBS})
