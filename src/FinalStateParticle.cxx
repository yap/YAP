#include "FinalStateParticle.h"

#include "Exceptions.h"
#include "logging.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
    : Particle(q, m, name),
      Model_(nullptr)
{
    // final state particles have fixed mass
    mass()->setVariableStatus(kFixed);
}

//-------------------------
bool FinalStateParticle::consistent() const
{
    bool C = Particle::consistent();

    if (ParticleCombinations_.empty()) {
        FLOG(ERROR) << "ParticleCombinations_ are empty!";
        C &= false;
    }

    return C;
}

//-------------------------
void FinalStateParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // pc must be final state particle
    if (!pc->isFinalStateParticle())
        throw exceptions::Exception("pc is not final state particle", "FinalStateParticle::addParticleCombination");

    // look for pc in ParticleCombinations_
    auto it = std::find_if(ParticleCombinations_.begin(), ParticleCombinations_.end(),
    [&](const std::shared_ptr<ParticleCombination>& p) {return p->indices() == pc->indices();});
    // if pc already contained, do nothing
    if (it == ParticleCombinations_.end())
        ParticleCombinations_.push_back(pc);
}

}

