#include "DecayChannel.h"
#include "Resonance.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude)
    :   Daughters_( {daughterA, daughterB}),
L_(L),
BlattWeisskopf_(this),
SpinAmplitude_(spinAmplitude),
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
    bool consistent = true;

    consistent &= this->daughterA()->consistent();
    consistent &= this->daughterB()->consistent();
    consistent &= BlattWeisskopf_.consistent();
    consistent &= SpinAmplitude_.consistent();


    // check if QuantumNumbers of SpinAmplitude objects match with Particles
    if (this->spinAmplitude().finalQuantumNumbersA() != this->daughterA()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of daughterA and SpinResonance don't match.";
        consistent = false;
    }

    if (this->spinAmplitude().finalQuantumNumbersB() != this->daughterB()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of daughterB and SpinResonance don't match.";
        consistent =  false;
    }

    // check if BlattWeisskopf points back to this DecayChannel
    if (this != BlattWeisskopf_.decayChannel()) {
        LOG(ERROR) << "DecayChannel::consistent() - BlattWeisskopf does not point back to this DecayChannel.";
        consistent =  false;
    }


    // check angular momentum conservation laws
    unsigned char l = this->l();
    unsigned char L_A = this->daughterA()->quantumNumbers().J();
    unsigned char L_B = this->daughterB()->quantumNumbers().J();

    if (l < abs(L_A - L_B) || l > abs(L_A + L_B)) {
        LOG(ERROR) << "DecayChannel::consistent() - angular momentum conservation violated.";
        consistent =  false;
    }

    // check if INITIAL QuantumNumbers of SpinAmplitude objects match with this Resonance's QuantumNumbers
    if (this->spinAmplitude().initialQuantumNumbers() != this->resonance()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of SpinAmplitude  and this channel's resonance don't match.";
        consistent =  false;
    }

    // check masses
    if (this->daughterA()->mass() + this->daughterB()->mass() > this->resonance()->mass()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        consistent =  false;
    }

    return consistent;
}

}
