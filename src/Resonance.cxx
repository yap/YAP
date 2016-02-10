#include "Resonance.h"

#include "DecayChannel.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "MassShape.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::shared_ptr<MassShape> massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(massShape)
{
    if (!MassShape_)
        throw exceptions::Exception("MassShape unset", "Resonance::Resonance");

    MassShape_->setResonance(this);
}

//-------------------------
std::complex<double> Resonance::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, unsigned dataPartitionIndex) const
{
    return DecayingParticle::amplitude(d, pc, two_m, dataPartitionIndex) * MassShape_->amplitude(d, pc, dataPartitionIndex);
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
void Resonance::addToModel()
{
    DecayingParticle::addToModel();
    MassShape_->addToModel();
}

//-------------------------
unsigned Resonance::addParticleCombination(std::shared_ptr<ParticleCombination> c)
{
    unsigned index = DecayingParticle::addParticleCombination(c);
    MassShape_->addParticleCombination(c);
    return index;
}

}
