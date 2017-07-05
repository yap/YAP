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

include_directories(
  ${YAP_SOURCE_DIR}/include
  ${source_dir}/include
  ${YAP_SOURCE_DIR}/tools
  ${YAP_SOURCE_DIR}/data
  )

set(YAP_TEST_SOURCES
  test_AmplitudeBasis.cxx
  test_Attributes.cxx
  test_BasisTransformations.cxx
  test_CalculationStatus.cxx
  test_ClebschGordan.cxx
  test_CompensatedSum.cxx
  test_DataSet.cxx
  test_deduce_parities.cxx
  test_FourMomenta.cxx
  test_FourMomentaCalculation.cxx
  test_Filter.cxx
  test_FreeAmplitude.cxx
  test_helicityFrame.cxx
  test_HelicityAngles.cxx
  test_HelicityAngles_boostRotate.cxx
  test_integration.cxx
  test_MathUtilities.cxx
  test_Matrix.cxx
  test_Model.cxx
  test_Parameter.cxx
  test_ParticleCombination.cxx
  test_ParticleTable.cxx
  test_Spin.cxx
  test_swapDalitzAxes.cxx
  test_swapFinalStates.cxx
  test_swapFourMomenta.cxx
  test_VariableStatus.cxx
  test_Vector.cxx
  test_WignerD.cxx
  yap_test.cxx
)

file(MAKE_DIRECTORY data)
configure_file(${YAP_SOURCE_DIR}/data/evt.pdl ${CMAKE_BINARY_DIR}/data/evt.pdl COPYONLY)

add_executable(test EXCLUDE_FROM_ALL ${YAP_TEST_SOURCES})
add_dependencies(test catch)
target_link_libraries(test YAP)
