#include "SpinAmplitude.h"

#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2)
    : InitialQuantumNumbers_(initial),
      FinalQuantumNumbers_( {{final1, final2}})
{}

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

    return consistent;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    return (InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_
            && FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]
            && FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]);
}

}
