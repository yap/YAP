#include "InitialStateParticle.h"

#include "CanonicalSpinAmplitude.h"
#include "logging.h"

#include <TLorentzRotation.h>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(q, mass, name, radialSize)
{
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    for (std::shared_ptr<yap::ParticleCombination>& pc : particleCombinations()) {
        pc->setParents();
    }
}


}
