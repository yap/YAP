#include "AmplitudeComponent.h"

#include "Parameter.h"

namespace yap {

//-------------------------
const VariableStatus RecalculableAmplitudeComponent::variableStatus() const
{
    return variable_status(parameters().begin(), parameters().end());
}

}
