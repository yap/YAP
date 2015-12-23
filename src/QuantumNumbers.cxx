#include "QuantumNumbers.h"

#include "logging.h"

namespace yap {

//-------------------------
bool QuantumNumbers::consistent() const
{
    bool result = true;

    // check charge parity for charged particle
    if (Q_ != 0 && C_ != 0) {
        FLOG(ERROR) << "charged particle has nonzero charge parity.";
        result = false;
    }

    /// \todo enable
    // check parity is set
    // if (P_ == 0) {
    //     LOG(ERROR) << "QuantumNumbers::consistent() - parity unset.";
    //     result = false;
    // }

    /*if (abs(TwoLambda_) > TwoJ_) {
        LOG(ERROR) << "QuantumNumbers::consistent() - Helicity is too big.";
        result = false;
    }*/

    return result;
}

//-------------------------
bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{
    return (lhs.twoJ() == rhs.twoJ() and lhs.P() == rhs.P() and lhs.C() == rhs.C() and
            lhs.twoI() == rhs.twoI() and lhs.G() == rhs.G() and lhs.Q() == rhs.Q());
}

}
