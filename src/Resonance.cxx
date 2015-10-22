#include "Resonance.h"

#include "DecayChannel.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::unique_ptr<MassShape>& massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(nullptr)
{
    setMassShape(massShape);
    Amplitude_->addDependencies(MassShape_->ParametersItDependsOn());
    Amplitude_->addDependencies(MassShape_->CachedDataValuesItDependsOn());
}

//-------------------------
void Resonance::setMassShape(std::unique_ptr<MassShape>& massShape)
{
    MassShape_.swap(massShape);
    MassShape_->borrowParametersFromResonance(this);
}

//-------------------------
std::complex<double> Resonance::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo check
    unsigned symIndex = symmetrizationIndex(pc);

    if (calculationStatus(pc, symIndex, d.index()) == kUncalculated) {

        std::complex<double> a = Complex_0;

        for (auto& c : channels()) {
            if (c->hasSymmetrizationIndex(pc))
                a += c->amplitude(d, pc);
        }

        a *= MassShape_->amplitude(d, pc);

        Amplitude_->setValue(a, d.dataPoint(), symIndex, d.index());

        DEBUG("Resonance " << name() << ": amplitude for " << std::string(*pc) << " = " << a);
    }

    return Amplitude_->value(d.dataPoint(), symIndex);
}

//-------------------------
bool Resonance::consistent() const
{
    bool consistent = true;

    consistent &= DecayingParticle::consistent();
    consistent &= MassShape_->consistent();

    if (! consistent) {
        LOG(ERROR) << "Resonance is not consistent:  " << this->name() << "\n";
    }

    return consistent;
}

//-------------------------
void Resonance::setInitialStateParticle(InitialStateParticle* isp)
{
    DecayingParticle::setInitialStateParticle(isp);
    // hand ISP to mass shape
    MassShape_->setInitialStateParticle(initialStateParticle());
}

//-------------------------
void Resonance::addChannel(std::unique_ptr<DecayChannel>& c)
{
    for (auto& pc : c->particleCombinations())
        MassShape_->addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(pc));

    DecayingParticle::addChannel(c);
}

//-------------------------
void Resonance::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    DecayingParticle::addSymmetrizationIndex(c);
    MassShape_->addSymmetrizationIndex(c);
}

//-------------------------
void Resonance::clearSymmetrizationIndices()
{
    DecayingParticle::clearSymmetrizationIndices();
    MassShape_->clearSymmetrizationIndices();
}


}
