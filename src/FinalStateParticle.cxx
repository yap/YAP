#include "FinalStateParticle.h"

#include "logging.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
    : Particle(q, m, name),
      InitialStateParticle_(nullptr)
{
    // final state particles have fixed mass
    mass()->setVariableStatus(kFixed);
}

//-------------------------
bool FinalStateParticle::consistent() const
{
    bool C = Particle::consistent();

    if (SymmetrizationIndices_.empty()) {
        FLOG(ERROR) << "SymmetrizationIndices_ are empty!";
        C &= false;
    }

    for (auto& pc : SymmetrizationIndices_)
        if (pc->indices().size() != 1) {
            FLOG(ERROR) << "ParticleCombination doesn't have size 1!";
            C &= false;
        }

    return C;
}

//-------------------------
void FinalStateParticle::setSymmetrizationIndexParents()
{
    if (!initialStateParticle())
        throw exceptions::InitialStateParticleUnset();

    ParticleCombinationVector PCs = SymmetrizationIndices_;

    // check if already set
    if (PCs[0]->parent())
        return;

    SymmetrizationIndices_.clear();

    for (auto& PC : PCs)
        for (auto& pc : initialStateParticle()->particleCombinationCache)
            if (ParticleCombination::equivDown(PC, pc.lock()))
                SymmetrizationIndices_.push_back(pc.lock());
}

}

