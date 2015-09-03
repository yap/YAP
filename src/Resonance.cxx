#include "Resonance.h"

#include "DecayChannel.h"
#include "FinalStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, MassShape* massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(massShape)
{}

//-------------------------
Amp Resonance::amplitude(DataPoint& d)
{
    // \todo implement
    return Amp(0.);
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
void Resonance::addChannel(DecayChannel* c)
{
    DecayingParticle::addChannel(c);

    for (std::shared_ptr<ParticleCombination> pc : c->particleCombinations()) {
        MassShape_->addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(pc));
        //MassShape_->addSymmetrizationIndex({pc});
    }
}

}
