#include "CachedValue.h"

#include "DataPoint.h"
#include "FourMomenta.h"
#include "logging.h"

namespace yap {

//-------------------------
CachedValue::CachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn,
                         std::vector<std::shared_ptr<CachedValue> > CachedValuesItDependsOn) :
    CalculationStatus_(kUncalculated)
{
    addDependencies(ParametersItDependsOn);
    addDependencies(CachedValuesItDependsOn);
}

//-------------------------
void CachedValue::setValue(std::complex<double> val)
{
    if (val != CachedValue_) {
        CachedValue_ = val;
        CalculationStatus_ = kCalculated;
    }
}

//-------------------------
void CachedValue::addDependencies(std::vector<std::shared_ptr<Parameter> > dep)
{
    for (auto& v : dep)
        ParametersItDependsOn_.insert(v);
}

//-------------------------
void CachedValue::addDependencies(std::vector<std::shared_ptr<CachedValue> > dep)
{
    for (auto& v : dep)
        CachedValuesItDependsOn_.insert(v);
}

//-------------------------
CalculationStatus CachedValue::calculationStatus()
{
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged)
            CalculationStatus_ = kUncalculated;
    }

    for (auto& v : CachedValuesItDependsOn_) {
        if (v->calculationStatus() == kUncalculated)
            CalculationStatus_ = kUncalculated;
    }

    return CalculationStatus_;
}

//-------------------------
void CachedValue::finishedPrecalculation()
{
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged)
            p->setVariableStatus(kUnchanged);
    }
}

}
