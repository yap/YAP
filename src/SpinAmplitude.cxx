#include "SpinAmplitude.h"

#include "logging.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial,
                             const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL)
    : AmplitudeComponentDataAccessor(isp),
      InitialQuantumNumbers_(initial),
      FinalQuantumNumbers_( {{final1, final2}}),
TwoL_(twoL)
{
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
    int twoJ_P = InitialQuantumNumbers_.twoJ();
    int twoJ_A = FinalQuantumNumbers_[0].twoJ();
    int twoJ_B = FinalQuantumNumbers_[1].twoJ();

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int twoL_AB = abs(twoJ_A - twoJ_B); twoL_AB <= abs(twoJ_A + twoJ_B); twoL_AB += 2) {
        for (int rhs = abs(TwoL_ - twoL_AB); rhs <= abs(TwoL_ + twoL_AB); rhs += 2) {
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
                   "J(parent) = " << spinToString(twoJ_P) << "; J(daughter1) = " << spinToString(twoJ_A) << "; J(daughter2) = " << spinToString(twoJ_B) << "; l = " << spinToString(TwoL_);
        consistent =  false;
    }

    return consistent;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    return (SymmetrizationIndices_ == rhs.SymmetrizationIndices_
            and InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_
            and FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]
            and FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]
            and TwoL_ == rhs.TwoL_);
}

}
