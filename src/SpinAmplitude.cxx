#include "SpinAmplitude.h"

#include "Constants.h"
#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char l)
    : InitialQuantumNumbers_(initial),
      FinalQuantumNumbers_( {{final1, final2}}),
L_(l)
{}

//-------------------------
Amp SpinAmplitude::amplitude(DataPoint& d)
{
    // \todo implement
    return Complex_0;
}

//-------------------------
bool SpinAmplitude::consistent() const
{
    bool consistent = true;

    // check charge conservation
    if (InitialQuantumNumbers_.Q() != FinalQuantumNumbers_[0].Q() + FinalQuantumNumbers_[1].Q()) {
        LOG(ERROR) << "SpinAmplitude::consistent() - charge conservation violated. " <<
                   "Q(parent) = " << (int)InitialQuantumNumbers_.Q() <<
                   "; Q(daughterA) = " << (int)FinalQuantumNumbers_[0].Q() <<
                   "; Q(daughterB) = " << (int)FinalQuantumNumbers_[1].Q();
        consistent =  false;
    }

    // check angular momentum conservation laws
    int twoL = 2 * L_;
    int twoJ_P = InitialQuantumNumbers_.twoJ();
    int twoJ_A = FinalQuantumNumbers_[0].twoJ();
    int twoJ_B = FinalQuantumNumbers_[1].twoJ();

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
        LOG(ERROR) << "SpinAmplitude::consistent() - angular momentum conservation violated. " <<
                   "J(parent) = " << .5 * twoJ_P << "; J(daughterA) = " << .5 * twoJ_A << "; J(daughterB) = " << .5 * twoJ_B << "; l = " << .5 * twoL;
        consistent =  false;
    }

    return consistent;
}

//-------------------------
bool operator== (const SpinAmplitude& lhs, const SpinAmplitude& rhs)
{
    return (lhs.InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_
            && lhs.FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]
            && lhs.FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]
            && lhs.L_ == rhs.L_);
}

}
