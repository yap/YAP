#include "CachedValue.h"

#include "logging.h"

namespace yap {

//-------------------------
CachedValueBase::CachedValueBase(ParameterSet pars) :
    ParametersItDependsOn_(pars),
    CalculationStatus_(kUncalculated)
{
}

//-------------------------
void CachedValueBase::removeDependency(std::shared_ptr<ParameterBase> dep)
{
    // look for parameter in set of dependencies
    auto it = ParametersItDependsOn_.find(dep);
    // if found, erase
    if (it != ParametersItDependsOn_.end())
        ParametersItDependsOn_.erase(it);
}
a
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
