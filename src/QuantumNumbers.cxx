#include "QuantumNumbers.h"

#include "logging.h"

namespace yap {

//-------------------------
QuantumNumbers::QuantumNumbers(unsigned char twoJ, char P, char C, unsigned char twoI, char G, char Q) :
    TwoJ_(twoJ),
    P_(P),
    C_(C),
    TwoI_(twoI),
    G_(G),
    Q_(Q)
{
}

//-------------------------
bool QuantumNumbers::consistent() const
{
    bool result = true;

    // check charge parity for charged particle
    if (Q_ != 0 && C_ != 0) {
        LOG(ERROR) << "QuantumNumbers::consistent() - charged particle has nonzero charge parity.";
        result = false;
    }

    // check parity is set
    if (P_ == 0) {
        LOG(ERROR) << "QuantumNumbers::consistent() - parity unset.";
        result = false;
    }

    return result;
}

//-------------------------
bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{
    //std::cout << lhs << " == " << rhs << "?\n";
    return (lhs.TwoJ_ == rhs.TwoJ_
            && lhs.P_ == rhs.P_
            && lhs.C_ == rhs.C_
            && lhs.TwoI_ == rhs.TwoI_
            && lhs.G_ == rhs.G_
            && lhs.Q_ == rhs.Q_);
}

//-------------------------
std::ostream& operator<< (std::ostream& os, const QuantumNumbers& obj)
{
    if (obj.C() != 0)
        os << "JPC = " << (int)obj.J() << (obj.P() > 0 ? "+" : "-") << (obj.C() > 0 ? "+" : "-");
    else
        os << "JP = " << (int)obj.J() << (obj.P() > 0 ? "+" : "-");

    //os << "\n";
    return os;
}

}
