set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
  BlattWeisskopf.cxx
  BreitWigner.cxx
  Constants.cxx
  DataAccessor.cxx
  DataSet.cxx
  DecayChannel.cxx
  DecayingParticle.cxx
  InitialStateParticle.cxx
  MassShape.cxx
  ParticleCombination.cxx
  Particle.cxx
  ParticleCombination.cxx
  ParticleFactory.cxx
  QuantumNumbers.cxx
  Resonance.cxx
  SpinAmplitude.cxx
)

make_static_library(YAP "${YAP_SOURCES}" "${ROOT_LIBS}")
#make_shared_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
