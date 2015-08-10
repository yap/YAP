#include "Resonance.h"
#include "FinalStateParticle.h"
#include "logging.h"

namespace yap {

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
