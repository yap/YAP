#include "CachedValueBase.h"

#include "logging.h"

namespace yap {

//-------------------------
CachedValueBase::CachedValueBase(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn) :
    CalculationStatus_(kUncalculated)
{
    addDependencies(ParametersItDependsOn);
}

//-------------------------
CalculationStatus CachedValueBase::calculationStatus()
{
    // if uncalculated, return without further checking
    if (CalculationStatus_ == kUncalculated)
        return CalculationStatus_;

    // else check if any dependencies are changed
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged) {
            // if so, update to uncalculated and return
            CalculationStatus_ = kUncalculated;
            return CalculationStatus_;
        }
    }

    // else return (calculated)
    return CalculationStatus_;
}

//-------------------------
void CachedValueBase::finishedPrecalculation()
{
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged)
            p->setVariableStatus(kUnchanged);
    }
}


}
