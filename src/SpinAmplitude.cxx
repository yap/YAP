#include "SpinAmplitude.h"

#include "Constants.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2)
    : InitialQuantumNumbers_(initial),
      FinalQuantumNumbers_( {{final1, final2}})
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
  // \todo implement
  return true;
}

}
