#include "InitialStateParticle.h"

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(q, mass, name, radialSize),
    DataAccessor()
{
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    bool result = true;

    result &= DecayingParticle::consistent();

    return result;
}

}
