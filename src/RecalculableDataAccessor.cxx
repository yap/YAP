#include "RecalculableDataAccessor.h"

#include "Parameter.h"

namespace yap {

//-------------------------
void RecalculableDataAccessor::setParameterFlagsToUnchanged()
{
    for (auto& p : parameters_)
        if (p->variableStatus() == VariableStatus::changed)
            p->setVariableStatus(VariableStatus::unchanged);
}

}

