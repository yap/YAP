#include "MassShapeWithNominalMass.h"

#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"

namespace yap {

//-------------------------
CalculationStatus MassShapeWithNominalMass::updateCalculationStatus(DataPartition& D) const
{
    // check if Resonance's mass has changed
    if (mass()->variableStatus() == VariableStatus::changed) {
        // if so, set calculationStatus to uncalculated for every particleCombination
        for (const auto& pc_symIndex : symmetrizationIndices()) {
            D.status(*T(), pc_symIndex.second) = CalculationStatus::uncalculated;
        }
        return CalculationStatus::uncalculated;
    }

    return CalculationStatus::calculated;
}

//-------------------------
void MassShapeWithNominalMass::setParameterFlagsToUnchanged()
{
    if (mass()->variableStatus() == VariableStatus::changed)
        mass()->setVariableStatus(VariableStatus::unchanged);
}

//-------------------------
std::shared_ptr<RealParameter> MassShapeWithNominalMass::mass()
{
    if (!resonance())
        throw exceptions::ResonanceUnset("MassShapeWithNominalMass::mass");
    return resonance()->mass();
}

//-------------------------
void MassShapeWithNominalMass::setParameters(const ParticleTableEntry& entry)
{
    try {
        mass()->setValue(entry.Mass);
    } catch (const exceptions::ResonanceUnset&) { /* ignore */ }
}

}




