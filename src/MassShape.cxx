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
CalculationStatus MassShape::updateCalculationStatus(DataPartition& d, std::shared_ptr<const ParticleCombination> c) const
{
    /// \todo implement; make clever
    /*for (ParameterStatus stat : parameterStatuses()) {
        if (stat == kChanged) {
            setCalculationStatus(c, kUncalculated);
            return kUncalculated;
        }
    }*/

    return kCalculated;
}

}
