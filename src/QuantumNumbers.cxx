#include "QuantumNumbers.h"

namespace yap {

//-------------------------
QuantumNumbers::QuantumNumbers(unsigned char twoJ, char P, char C, char I, char G) :
    TwoJ_(twoJ),
    P_(P),
    C_(C),
    I_(I),
    G_(G)
{}

//-------------------------
bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{
    return (lhs.TwoJ_ == rhs.TwoJ_
            && lhs.P_ == rhs.P_
            && lhs.C_ == rhs.C_
            && lhs.I_ == rhs.I_
            && lhs.G_ == rhs.G_);
}

//-------------------------
std::ostream& operator<< (std::ostream& os, const QuantumNumbers& obj) {
       os << "J: " << obj.J() << " P: " << obj.P() << " C: " << obj.C() << "\n";
       return os;
}

}
