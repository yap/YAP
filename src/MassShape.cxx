#include "MassShape.h"

#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape(InitialStateParticle* isp) :
    AmplitudeComponentDataAccessor(isp, &ParticleCombination::equivByOrderlessContent),
    ParameterSet()
{
}

//-------------------------
CalculationStatus MassShape::updateCalculationStatus(std::shared_ptr<ParticleCombination> c)
{
    for (ParameterStatus stat : parameterStatuses()) {
        if (stat == kChanged) {
            setCalculationStatus(c, kUncalculated);
            return kUncalculated;
        }
    }

    return kCalculated;
}

}
