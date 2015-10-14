#include "CachedDataValue.h"

#include "logging.h"

namespace yap {

//-------------------------
CachedDataValueBase::CachedDataValueBase(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn,
        std::vector<std::shared_ptr<CachedDataValueBase> > CachedDataValuesItDependsOn)
{
    addDependencies(ParametersItDependsOn);
    addDependencies(CachedDataValuesItDependsOn);
}

//-------------------------
CalculationStatus CachedDataValueBase::calculationStatus(unsigned symmetrizationIndex, unsigned dataPartitionIndex)
{
    // if uncalculated, return without further checking
    if (CalculationStatus_[dataPartitionIndex][symmetrizationIndex] == kUncalculated)
        return kUncalculated;

    // else check if any dependencies are changed
    for (auto& c : CachedDataValuesItDependsOn_) {
        if (c->variableStatus(symmetrizationIndex, dataPartitionIndex) == kChanged) {
            CalculationStatus_[dataPartitionIndex][symmetrizationIndex] = kUncalculated;
            return kUncalculated;
        }
    }
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged) {
            // if so, update to uncalculated and return
            CalculationStatus_[dataPartitionIndex][symmetrizationIndex] = kUncalculated;
            return kUncalculated;
        }
    }

    // else return (calculated)
    return CalculationStatus_[dataPartitionIndex][symmetrizationIndex];
}

}
