#include "DecayChannel.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle& daughterA, Particle& daughterB, unsigned int L, SpinAmplitude& spinAmplitude)
:   L_(L), spinAmplitude_(spinAmplitude),
    freeAmplitude_(0)
{
  std::get<0>(daughters_) = daughterA;
  std::get<1>(daughters_) = daughterB;
}

//-------------------------
Amp DecayChannel::amplitude(DataPoint& d) {
  return Amp(1);
}

//-------------------------
bool DecayChannel::checkConsistency() const {
  return true;
}

//-------------------------
const Particle& getDaughter(unsigned int i) {
  switch (i) {
  case 0:
    return std::get<0>(daughters_);
  case 1:
    return std::get<1>(daughters_);
  default:

  }

}

}
