#include "FinalStateParticle.h"

#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
    : Particle(q, m, name),
      Model_(nullptr)
{
    // final state particles have fixed mass
    mass()->variableStatus() = VariableStatus::fixed;
}

//-------------------------
void FinalStateParticle::addParticleCombination(const std::shared_ptr<ParticleCombination>& pc)
{
    // pc must be final state particle
    if (!pc->isFinalStateParticle())
        throw exceptions::Exception("pc is not final state particle", "FinalStateParticle::addParticleCombination");

    // if pc not already in particleCombinations vector, add it
    if (!any_of(particleCombinations(), pc))
        Particle::addParticleCombination(pc);
}

//-------------------------
bool valid_final_state(const std::shared_ptr<ParticleCombination>& pc, const FinalStateParticleVector& FSPs)
{
    if (pc->indices().size() != FSPs.size())
        return false;

    auto I = pc->indices();
    for (const auto& fsp : FSPs) {
        auto it = I.end();
        for (const auto& pc : fsp->particleCombinations()) {
            it = std::find(I.begin(), I.end(), pc->indices()[0]);
            if (it != I.end())
                break;
        }
        if (it == I.end())
            return false;
        I.erase(it);
    }
    return true;
}

//-------------------------
bool is_final_state_particle(const std::shared_ptr<Particle>& p)
{
    return std::dynamic_pointer_cast<FinalStateParticle>(p) != nullptr;
}


}

