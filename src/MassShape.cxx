#include "MassShape.h"

#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape(std::initializer_list<double> pars) :
    AmplitudeComponentDataAccessor(&ParticleCombination::equivByOrderlessContent),
    ParameterSet(pars)
{
}

//-------------------------
bool MassShape::setParameters(const ParticleTableEntry& entry)
{
    if (entry.MassShapeParameters_.size() < Parameters_.size()) {
        LOG(ERROR) << "MassShape::setParameters - entry.MassShapeParameters_ is too small for MassShape object.";
        return false;
    }
    Parameters_.assign(entry.MassShapeParameters_.begin(), entry.MassShapeParameters_.begin() + Parameters_.size());
    return true;
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
