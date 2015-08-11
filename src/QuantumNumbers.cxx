#include "QuantumNumbers.h"

namespace yap {

//-------------------------
QuantumNumbers::QuantumNumbers(unsigned char twoJ, signed char P, signed char C, signed char I, signed char G) :
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
  if (obj.C() != 0)
    os << "JPC = " << (int)obj.J() << (obj.P()>0?"+":"-") << (obj.C()>0?"+":"-");
  else
    os << "JP = " << (int)obj.J() << (obj.P()>0?"+":"-");

  os << "\n";
  return os;
}

}
