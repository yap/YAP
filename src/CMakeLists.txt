set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
  DataPoint.cxx
)

make_static_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
#make_shared_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
