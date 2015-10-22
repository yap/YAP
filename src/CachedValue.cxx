#include "CachedValue.h"

#include "logging.h"

namespace yap {

//-------------------------
CachedValueBase::CachedValueBase(std::vector<std::shared_ptr<ComplexParameter> > ParametersItDependsOn) :
    CalculationStatus_(kUncalculated)
{
    addDependencies(ParametersItDependsOn);
}

//-------------------------
void CachedValueBase::removeDependency(std::shared_ptr<ComplexParameter> dep)
{
    // look for parameter in set of dependencies
    auto it = ParametersItDependsOn_.find(dep);
    // if found, erase
    if (it != ParametersItDependsOn_.end())
        ParametersItDependsOn_.erase(it);
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
