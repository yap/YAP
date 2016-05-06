#include "CachedValue.h"

#include "Parameter.h"

namespace yap {

//-------------------------
CachedValueBase::CachedValueBase(ParameterSet pars) :
    ParametersItDependsOn_(pars),
    CalculationStatus_(CalculationStatus::uncalculated)
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

//-------------------------
const CalculationStatus& CachedValueBase::calculationStatus()
{
    // if uncalculated, return without further checking
    if (CalculationStatus_ == CalculationStatus::uncalculated)
        return CalculationStatus_;

    // else check if any dependencies are changed
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == VariableStatus::changed) {
            // if so, update to uncalculated and return
            CalculationStatus_ = CalculationStatus::uncalculated;
            return CalculationStatus_;
        }
    }

    // else return (calculated)
    return CalculationStatus_;
}

}
