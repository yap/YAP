#include "RecalculableDataAccessor.h"

#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void RecalculableDataAccessor::setParameterFlagsToUnchanged()
{
    for (auto& p : Parameters_)
        if (p->variableStatus() == VariableStatus::changed)
            p->variableStatus() = VariableStatus::unchanged;
}

//-------------------------
void RecalculableDataAccessor::registerWithModel()
{
    DataAccessor::registerWithModel();
    const_cast<Model*>(model())->RecalculableDataAccessors_.insert(this);
}

}

