#include "SpinAmplitude.h"

#include "logging.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial,
                             const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL)
    : AmplitudeComponentDataAccessor(),
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
    if (!angularMomentumConserved(InitialQuantumNumbers_, FinalQuantumNumbers_[0], FinalQuantumNumbers_[1], TwoL_)) {
        LOG(ERROR) << "SpinAmplitude::consistent() - angular momentum conservation violated. "
                   << "J(parent) = " << spinToString(InitialQuantumNumbers_.twoJ())
                   << "; J(daughter1) = " << spinToString(FinalQuantumNumbers_[0].twoJ())
                   << "; J(daughter2) = " << spinToString(FinalQuantumNumbers_[1].twoJ())
                   << "; l = " << spinToString(TwoL_);
        consistent =  false;
    }

    return consistent;
}

//-------------------------
bool SpinAmplitude::angularMomentumConserved(const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2,
        unsigned char twoL)
{
    int twoJ_P = initial.twoJ();
    int twoJ_A = final1.twoJ();
    int twoJ_B = final2.twoJ();

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

    return ok;
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
