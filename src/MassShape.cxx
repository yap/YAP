#include "MassShape.h"

#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape(InitialStateParticle* isp) :
    AmplitudeComponent(),
    DataAccessor(isp, &ParticleCombination::equivByOrderlessContent),
    ParameterSet()
{
}

}
