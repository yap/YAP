#include "Resonance.h"
#include "FinalStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, const MassShape& massShape) :
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
    consistent &= MassShape_.consistent();

    return consistent;
}

}
