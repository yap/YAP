#include "MassShape.h"

#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape(InitialStateParticle* isp) :
    AmplitudeComponentDataAccessor(isp, &ParticleCombination::equivByOrderlessContent),
    ParameterSet()
{
}

}
