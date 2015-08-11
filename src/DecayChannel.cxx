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
    unsigned char l = decayAngularMomentum();
    unsigned char L_A = Daughters_[0]->quantumNumbers().J();
    unsigned char L_B = Daughters_[1]->quantumNumbers().J();

    if (l < abs(L_A - L_B) || l > abs(L_A + L_B)) {
        LOG(ERROR) << "DecayChannel::consistent() - angular momentum conservation violated.";
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
