#include "InitialStateParticle.h"

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

}
