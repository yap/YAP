#include "RecalculableDataAccessor.h"

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
const VariableStatus variableStatus(const RecalculableDataAccessor& rda)
{
    for (auto& p : rda.parameters())
        if (p->variableStatus() == VariableStatus::changed)
            return VariableStatus::changed;

    return VariableStatus::unchanged;
}

}

