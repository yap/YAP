#include "DecayChannel.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude)
:   Daughters_({daughterA, daughterB}),
    L_(L), SpinAmplitude_(spinAmplitude),
    FreeAmplitude_(0)
{
  ;
}

//-------------------------
Amp DecayChannel::amplitude(DataPoint& d) {
  return Amp(1);
}

//-------------------------
bool DecayChannel::checkConsistency() const {
  return true;
}

}
