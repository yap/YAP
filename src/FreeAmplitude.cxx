#include "FreeAmplitude.h"

#include "DecayChannel.h"
#include "spin.h"
#include "SpinAmplitude.h"

namespace yap {

//-------------------------
std::string to_string(const FreeAmplitude& fa)
{
    return to_string(*fa.decayChannel())
        + ", M = " + spin_to_string(fa.twoM())
        + ", " + to_string(*fa.spinAmplitude());
}

}
