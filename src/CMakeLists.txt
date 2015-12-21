set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedDataValue.cxx
	CachedValue.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPoint.cxx
	DataSet.cxx
	DecayChannel.cxx
	DecayingParticle.cxx
	FinalStateParticle.cxx
	FourMomenta.cxx
	HelicityAngles.cxx
	HelicitySpinAmplitude.cxx
	InitialStateParticle.cxx
	MeasuredBreakupMomenta.cxx
	ParticleCombination.cxx
  ParticleCombinationCache.cxx
	Particle.cxx
	ParticleFactory.cxx
	QuantumNumbers.cxx
	Resonance.cxx
	SpinAmplitude.cxx
	SpinUtilities.cxx
	WignerD.cxx
)

make_static_library(YAP "${YAP_SOURCES}" "${ROOT_LIBS}")
#make_shared_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
