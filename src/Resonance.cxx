#include "Resonance.h"

#include "DecayChannel.h"
#include "FinalStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::shared_ptr<MassShape> massShape) :
    DecayingParticle(massShape->initialStateParticle(), q, mass, name, radialSize),
    MassShape_(massShape)
{}

//-------------------------
Amp Resonance::calcAmplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo check
    Amp a = Complex_0;

    for (auto& c : channels()) {
        if (c->hasSymmetrizationIndex(pc))
            a += c->freeAmplitude() * c->amplitude(d, pc);
    }

    a *= MassShape_->amplitude(d, pc);

    DEBUG("Resonance " << name() << ": amplitude for " << std::string(*pc) << " = " << a);

    return a;
}

//-------------------------
bool Resonance::consistent() const
{
    bool consistent = true;

    consistent &= DecayingParticle::consistent();
    consistent &= massShape()->consistent();

    if (! consistent) {
        LOG(ERROR) << "Resonance is not consistent:  " << this->name() << "\n";
    }

    return consistent;
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
