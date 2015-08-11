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

    // check daughters
    for (Daughters::const_iterator d = Daughters_.begin(); d != Daughters_.end(); ++d)
        consistent &= (*d)->consistent();

    consistent &= BlattWeisskopf_.consistent();
    consistent &= SpinAmplitude_.consistent();

    // check size of spin amplitude quantum numbers and size of daughters
    if (SpinAmplitude_.finalQuantumNumbers().size() != Daughters_.size()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum mbumbers object and daughters object size mismatch";
        consistent = false;
    }

    // check if QuantumNumbers of SpinAmplitude objects match with Particles
    for (unsigned i = 0; i < Daughters_.size(); ++i) {
        if (SpinAmplitude_.finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
            LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of " << i << " and SpinResonance don't match.";
            consistent = false;
        }
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
    if (spinAmplitude().initialQuantumNumbers() != parent()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of SpinAmplitude  and this channel's resonance don't match.";
        consistent =  false;
    }

    // check masses
    double finalMass = 0;
    for (Daughters::const_iterator d = Daughters_.begin(); d != Daughters_.end(); ++d)
        finalMass += (*d)->mass();
    if (finalMass > parent()->mass()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        consistent =  false;
    }

    return consistent;
}

}
