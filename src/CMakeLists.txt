set(INCLUDE_DIRECTORIES
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	AmplitudeComponent.cxx
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedValue.cxx
	ClebschGordan.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPoint.cxx
	DataSet.cxx
	DecayChannel.cxx
	DecayTree.cxx
	DecayTreeVectorIntegral.cxx
	DecayingParticle.cxx
	Filters.cxx
	FinalStateParticle.cxx
	Flatte.cxx
	FourMomenta.cxx
	FreeAmplitude.cxx
	HelicityAngles.cxx
	HelicityFormalism.cxx
	ImportanceSampler.cxx
	Integrator.cxx
	MassRange.cxx
	MassShape.cxx
	MassShapeWithNominalMass.cxx
	MeasuredBreakupMomenta.cxx
	Model.cxx
	ModelIntegral.cxx
	PDL.cxx
	Parameter.cxx
	Particle.cxx
	ParticleCombination.cxx
	ParticleCombinationCache.cxx
	ParticleFactory.cxx
	PoleMass.cxx
	QuantumNumbers.cxx
	RecalculableDataAccessor.cxx
	RelativisticBreitWigner.cxx
	Resonance.cxx
	SpinAmplitude.cxx
	SpinAmplitudeCache.cxx
	StaticDataAccessor.cxx
	StatusManager.cxx
	WignerD.cxx
	ZemachFormalism.cxx
)

add_library(YAP SHARED ${YAP_SOURCES})

# Install rule for libYAP.so
install(TARGETS YAP LIBRARY
        DESTINATION lib/
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ             GROUP_EXECUTE
                    WORLD_READ             WORLD_EXECUTE)
