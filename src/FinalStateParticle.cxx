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
    mass()->setVariableStatus(VariableStatus::fixed);
}

//-------------------------
void FinalStateParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // pc must be final state particle
    if (!pc->isFinalStateParticle())
        throw exceptions::Exception("pc is not final state particle", "FinalStateParticle::addParticleCombination");

    // look for pc in ParticleCombinations_
    auto it = std::find_if(particleCombinations().begin(), particleCombinations().end(),
    [&](const std::shared_ptr<ParticleCombination>& p) {return ParticleCombination::equalUpAndDown(p, pc);});
    // if pc already contained, do nothing
    if (it == particleCombinations().end())
        Particle::addParticleCombination(pc);
}

}

