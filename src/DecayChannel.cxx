#include "DecayChannel.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude)
    :   Daughters_( {daughterA, daughterB}),
L_(L), SpinAmplitude_(spinAmplitude),
FreeAmplitude_(0)
{
    ;
}

//-------------------------
Amp DecayChannel::amplitude(DataPoint& d)
{
    return Amp(1);
}

//-------------------------
bool DecayChannel::consistent() const
{

    // check if QuantumNumbers of SpinAmplitude objects match with Particles
    if (this->spinAmplitude().finalQuantumNumbersA() != this->daughterA()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of daughterA and SpinResonance don't match.";
        return false;
    }

    if (this->spinAmplitude().finalQuantumNumbersB() != this->daughterB()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of daughterB and SpinResonance don't match.";
        return false;
    }


    return true;
}

}
