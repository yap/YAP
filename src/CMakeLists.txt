set(INCLUDE_DIRECTORIES
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	AmplitudePair.cxx
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedDataValue.cxx
	CachedValue.cxx
	ClebschGordan.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPoint.cxx
	DataSet.cxx
	DecayChannel.cxx
	DecayingParticle.cxx
	FinalStateParticle.cxx
	Flatte.cxx
	FourMomenta.cxx
	HelicityAngles.cxx
	HelicityFormalism.cxx
	MassShape.cxx
	MassShapeWithNominalMass.cxx
	MeasuredBreakupMomenta.cxx
	Model.cxx
	ParticleCombination.cxx
	ParticleCombinationCache.cxx
	Particle.cxx
	ParticleFactory.cxx
	PoleMass.cxx
	QuantumNumbers.cxx
	Resonance.cxx
	SpinAmplitude.cxx
	StatusManager.cxx
	WignerD.cxx
	ZemachFormalism.cxx
)

add_library(YAP SHARED ${YAP_SOURCES})

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=<lib_path> <path_to_CMakeLists.tex>
# otherwise, default LD_LIBRARY_PATH
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set(LIBRARY_OUTPUT_DIRECTORY ${YAP_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE})
endif()

install(TARGETS YAP LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
