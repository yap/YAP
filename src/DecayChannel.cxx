#include "DecayChannel.h"
#include "Resonance.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L,
                           std::shared_ptr<SpinAmplitude> spinAmplitude) :
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
    if (! consistent) {
        LOG(ERROR) << "DecayChannel::consistent() - daughter inconsistent";
    }

    consistent &= BlattWeisskopf_.consistent();
    consistent &= spinAmplitude()->consistent();

    // check size of spin amplitude quantum numbers and size of daughters
    if (spinAmplitude()->finalQuantumNumbers().size() != Daughters_.size()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers object and daughters object size mismatch";
        consistent = false;
    }

    // check if QuantumNumbers of SpinAmplitude objects match with Particles
    if (spinAmplitude()->initialQuantumNumbers() != parent()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of parent and SpinResonance don't match.";
        consistent = false;
    }

    for (unsigned i = 0; i < Daughters_.size(); ++i) {
        if (spinAmplitude()->finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
            LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of " << i << " and SpinResonance don't match.";
            consistent = false;
        }
    }

    // check if BlattWeisskopf points back to this DecayChannel
    if (this != BlattWeisskopf_.decayChannel()) {
        LOG(ERROR) << "DecayChannel::consistent() - BlattWeisskopf does not point back to this DecayChannel.";
        consistent =  false;
    }

    // check charge conservation
    if (this->parent()->quantumNumbers().Q() != this->daughters()[0]->quantumNumbers().Q() + this->daughters()[1]->quantumNumbers().Q()) {
        LOG(ERROR) << "DecayChannel::consistent() - charge conservation violated. " <<
                   "Q(parent) = " << (int)this->parent()->quantumNumbers().Q() <<
                   "; Q(daughterA) = " << (int)this->daughters()[0]->quantumNumbers().Q() <<
                   "; Q(daughterB) = " << (int)this->daughters()[1]->quantumNumbers().Q();
        consistent =  false;
    }

    // check angular momentum conservation laws
    int twoL = 2 * this->decayAngularMomentum();
    int twoJ_P = this->parent()->quantumNumbers().twoJ();
    int twoJ_A = this->daughters()[0]->quantumNumbers().twoJ();
    int twoJ_B = this->daughters()[1]->quantumNumbers().twoJ();

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int twoL_AB = abs(twoJ_A - twoJ_B); twoL_AB <= abs(twoJ_A + twoJ_B); twoL_AB += 2) {
        for (int rhs = abs(twoL - twoL_AB); rhs <= abs(twoL + twoL_AB); rhs += 2) {
            if (twoJ_P == rhs) {
                ok = true;
                break;
            }
            if (ok)
                break;
        }
    }

    if (!ok) {
        LOG(ERROR) << "DecayChannel::consistent() - angular momentum conservation violated. " <<
                   "J(parent) = " << .5 * twoJ_P << "; J(daughterA) = " << .5 * twoJ_A << "; J(daughterB) = " << .5 * twoJ_B << "; l = " << .5 * twoL;
        consistent =  false;
    }

    // check if INITIAL QuantumNumbers of SpinAmplitude objects match with this Resonance's QuantumNumbers
    if (spinAmplitude()->initialQuantumNumbers() != parent()->quantumNumbers()) {
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

    if (! consistent) {
        LOG(ERROR) << "Channel is not consistent:  " << this->parent()->name() << " - > " << this->daughters()[0]->name() << " + " << this->daughters()[1]->name() << "\n";
    }
    return consistent;
}

}
