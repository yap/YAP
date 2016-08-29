#include "RecalculableDataAccessor.h"

#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void RecalculableDataAccessor::setParameterFlagsToUnchanged()
{
    for (auto& p : parameters_)
        if (p->variableStatus() == VariableStatus::changed)
            p->variableStatus() = VariableStatus::unchanged;
}

//-------------------------
const VariableStatus variable_status(const RecalculableDataAccessor& rda)
{
    for (auto& p : rda.parameters())
        if (p->variableStatus() == VariableStatus::changed)
            return VariableStatus::changed;

    return VariableStatus::unchanged;
}

//-------------------------
void RecalculableDataAccessor::registerWithModel()
{
    DataAccessor::registerWithModel();
    const_cast<Model*>(model())->RecalculableDataAccessors_.insert(this);
}

}

