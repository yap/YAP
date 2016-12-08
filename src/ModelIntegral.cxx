#include "ModelIntegral.h"

#include "DecayingParticle.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "Model.h"
#include "Parameter.h"

#include <memory>

namespace yap {

//-------------------------
ModelComponentIntegral::ModelComponentIntegral(const ModelComponent& c)
    : Admixture(c.admixture()),
      Integral(c.decayTrees())
{
}
    
//-------------------------
ModelIntegral::ModelIntegral(const Model& model)
{
    Integrals_.reserve(model.components().size());
    for (const auto& c : model.components())
        Integrals_.emplace_back(c);
}

//-------------------------
const RealIntegralElement integral(const ModelIntegral& MI)
{
    return std::accumulate(MI.integrals().begin(), MI.integrals().end(), RealIntegralElement(),
                           [](RealIntegralElement& I, const ModelComponentIntegral& mci)
                           { return I += mci.Admixture->value() * integral(mci.Integral); });
}

}
