#include "Resonance.h"
#include "FinalStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Amp amplitude(DataPoint& d)
{
    // TODO implement
    return Amp(1);
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
