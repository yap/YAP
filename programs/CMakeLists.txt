#include_directories(${BellePWA_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${BAT_INCLUDE_DIR})
include_directories(${YAP_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

#make_executable(test test.cxx BellePWA "${ROOT_LIBS}" "{BAT_LIBS}")
make_executable(D4piTest D4piTest YAP "${ROOT_LIBS}")
make_executable(D3piTest D3piTest YAP "${ROOT_LIBS}")
make_executable(DKKpiTest DKKpiTest YAP "${ROOT_LIBS}")

#make_executable(GeneratorTest GeneratorTest YAP "${ROOT_LIBS}")

