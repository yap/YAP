#include "ModelIntegral.h"

#include "DecayingParticle.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "Model.h"
#include "Parameter.h"

#include <memory>

namespace yap {

//-------------------------
ModelIntegral::ModelIntegral(const Model& model)
{
    // for each initial state particle
    for (const auto& isp_mix : model.initialStateParticles())
        // for each spin projection
        for (const auto& m_b : isp_mix.second)
            // create new DecayTreeVectorIntegral
            Integrals_.emplace(m_b.second, DecayTreeVectorIntegral(isp_mix.first->decayTrees().at(m_b.first)));
}

//-------------------------
const RealIntegralElement ModelIntegral::integral() const
{
   return std::accumulate(Integrals_.begin(), Integrals_.end(), RealIntegralElement(),
                          [](RealIntegralElement& i, const IntegralMap::value_type& b_I)
                          { return i += b_I.first->value() * b_I.second.integral(); });
}

//-------------------------
const DecayTreeVectorIntegral& ModelIntegral::integral(const DecayTreeVector& dtv) const
{
    auto it = find_if(Integrals_.begin(), Integrals_.end(),
                      [&](const IntegralMap::value_type& b_I)
                      {return b_I.second.decayTrees() == dtv;});
    if (it == Integrals_.end())
        throw exceptions::Exception("DecayTreeVector not found", "ModelIntegral::integral");
    return it->second;
}

//-------------------------
const Model* ModelIntegral::model() const
{
    return Integrals_.empty() ? nullptr : Integrals_.begin()->second.model();
}

}
