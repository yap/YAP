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
    setVariableStatus(VariableStatus::changed);
}

//-------------------------
void set_values(ParameterVector& pars, const std::vector<double>& vals)
{
    auto val_it = vals.begin();
    auto par_it = pars.begin();
    for (; val_it != vals.end() and par_it != pars.end();  ++val_it, ++par_it)
        val_it = set_value(**par_it, val_it);
    // if either iterator is not at end, then there was a size mismatch
    if (val_it != vals.end() or par_it != pars.end())
        throw exceptions::Exception("size mismatch", "set_values");
}


}
