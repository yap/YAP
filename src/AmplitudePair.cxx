#include "AmplitudePair.h"

#include "CachedDataValue.h"
#include "DecayChannel.h"

namespace yap {

//-------------------------
AmplitudePair::AmplitudePair(DecayChannel* dc, std::complex<double> free) :
    Fixed(ComplexCachedDataValue::create(dc)),
    Free(std::make_shared<ComplexParameter>(free))
{
}

}
