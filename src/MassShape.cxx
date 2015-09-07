#include "MassShape.h"

#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    AmplitudeComponent(),
    DataAccessor(ParticleCombination::equivByOrderlessContent),
    ParameterSet()
{
}

}
