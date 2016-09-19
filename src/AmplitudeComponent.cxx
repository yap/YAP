#include "AmplitudeComponent.h"

#include "Parameter.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
const bool StaticAmplitudeComponent::validFor(const ParticleCombination& pc) const
{
    return symmetrizationIndices().find(pc.shared_from_this()) != symmetrizationIndices().end();
}

//-------------------------
const bool RecalculableAmplitudeComponent::validFor(const ParticleCombination& pc) const
{
    return symmetrizationIndices().find(pc.shared_from_this()) != symmetrizationIndices().end();
}

//-------------------------
const VariableStatus RecalculableAmplitudeComponent::status() const
{
    return variable_status(parameters().begin(), parameters().end());
}

}
