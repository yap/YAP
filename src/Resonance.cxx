#include "Resonance.h"

#include "DecayChannel.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::unique_ptr<MassShape> massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(std::move(massShape))
{
    if (!MassShape_)
        throw exceptions::Exception("MassShape unset", "Resonance::Resonance");

    MassShape_->setResonance(this);
}

//-------------------------
bool Resonance::consistent() const
{
    bool C = DecayingParticle::consistent();

    if (!MassShape_->consistent()) {
        FLOG(ERROR) << "MassShape not consistent";
        C &= false;
    }
    if (MassShape_->resonance() != this) {
        FLOG(ERROR) << "MassShape's resonance is not this";
        C &= false;
    }

    return C;
}

//-------------------------
void Resonance::addChannel(std::unique_ptr<DecayChannel> c)
{
    DecayingParticle::addChannel(std::move(c));

    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "Resonance::addChannel");

    for (auto& pc : channels().back()->particleCombinations())
        MassShape_->addSymmetrizationIndex(pc);

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

//-------------------------
DataAccessorSet Resonance::dataAccessors()
{
    // call DecayingParticle's function
    DataAccessorSet V = DecayingParticle::dataAccessors();

    // add mass shape
    V.emplace(MassShape_);

    return V;
}

}
