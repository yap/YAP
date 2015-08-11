#include "DecayChannel.h"
#include "Resonance.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude) :
    Daughters_( {daughterA, daughterB}),
            L_(L),
            BlattWeisskopf_(this),
            SpinAmplitude_(spinAmplitude),
            FreeAmplitude_(0)
{}

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
    unsigned char s_P = this->parent()->quantumNumbers().J();
    unsigned char s_A = this->daughterA()->quantumNumbers().J();
    unsigned char s_B = this->daughterB()->quantumNumbers().J();

    std::cout << this->parent()->quantumNumbers();
    std::cout << "J(parent) = " << s_P << "; J(daughterA) = " << s_A << "; J(daughterB) = " << s_B << "; l = " << l << "\n";

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int l_AB = abs(s_A - s_B); l_AB <= abs(s_A + s_B); ++l_AB) {
      for (int rhs = abs(l - l_AB); rhs <= abs(l + l_AB); ++rhs) {
        if (s_P == rhs) {
          ok = true;
          break;
        }
        if (ok)
          break;
      }
    }

    if (!ok) {
        LOG(ERROR) << "DecayChannel::consistent() - angular momentum conservation violated. " <<
            "J(parent) = " << s_P << "; J(daughterA) = " << s_A << "; J(daughterB) = " << s_B << "; l = " << l;
        consistent =  false;
    }

    // check if INITIAL QuantumNumbers of SpinAmplitude objects match with this Resonance's QuantumNumbers
    if (this->spinAmplitude().initialQuantumNumbers() != this->parent()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of SpinAmplitude  and this channel's resonance don't match.";
        consistent =  false;
    }

    // check masses
    if (this->daughterA()->mass() + this->daughterB()->mass() > this->parent()->mass()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
void DecayChannel::setParent(DecayingParticle* parent)
{
    Parent_ = parent;
}

}
