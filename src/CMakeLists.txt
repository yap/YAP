set(INCLUDE_DIRECTORIES
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	AmplitudeComponent.cxx
	Attributes.cxx
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedValue.cxx
	ClebschGordan.cxx
	ConstantWidthBreitWigner.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPoint.cxx
	DataSet.cxx
	DecayChannel.cxx
	DecayingParticle.cxx
	DecayTree.cxx
	DecayTreeVectorIntegral.cxx
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
	NonrelativisticBreitWigner.cxx
	NonrelativisticConstantWidthBreitWigner.cxx
	Particle.cxx
	ParticleCombination.cxx
	ParticleCombinationCache.cxx
	ParticleTable.cxx
	PDL.cxx
	PoleMass.cxx
	RecalculableDataAccessor.cxx
	SpinAmplitude.cxx
	SpinAmplitudeCache.cxx
	StaticDataAccessor.cxx
	StatusManager.cxx
	UnitSpinAmplitude.cxx
	WignerD.cxx
	ZemachFormalism.cxx
)

add_library(YAP SHARED ${YAP_SOURCES})

# Install rule for libYAP.so
install(TARGETS YAP LIBRARY
        COMPONENT   core
        DESTINATION lib/
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ             GROUP_EXECUTE
                    WORLD_READ             WORLD_EXECUTE)
