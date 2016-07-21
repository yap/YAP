#include "Parameter.h"

namespace yap {

//-------------------------
ComplexComponentParameter::ComplexComponentParameter(std::shared_ptr<ComplexParameter> par) :
    RealParameter(),
    Parent_(par)
{
    if (!Parent_)
        throw exceptions::Exception("Parent unset", "RealSubparameter::RealSubparameter");
}


//-------------------------
void ComplexComponentParameter::setValue(double val)
{
    if (variableStatus() == VariableStatus::fixed)
        throw exceptions::ParameterIsFixed("", "ComplexComponentParameter::setValue");
    if (value() == val)
        return;
    Parent_->setValue(setComponent(Parent_->value(), val));
    variableStatus() = VariableStatus::changed;
}

}
