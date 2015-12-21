#include "FinalStateParticle.h"

#include "logging.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
    : Particle(q, m, name)
{
    // final state particles have fixed mass
    mass()->setVariableStatus(kFixed);
}

//-------------------------
bool FinalStateParticle::consistent() const
{
    bool consistent = true;

    consistent &= Particle::consistent();

    if (SymmetrizationIndices_.empty()) {
        LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ are empty!";
        return false;
    }

    for (auto i : SymmetrizationIndices_) {
        if (i->indices().size() != 1) {
            LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ don't have size 1!";
            return false;
        }
        if (i->daughters().size() != 0) {
            LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ have daughters!";
            return false;
        }
    }

    return consistent;
}

//-------------------------
void FinalStateParticle::setSymmetrizationIndexParents()
{
    ParticleCombinationVector PCs = SymmetrizationIndices_;

    // check if already set
    if (PCs[0]->parent())
        return;

    SymmetrizationIndices_.clear();

    for (auto& PC : PCs)
        for (auto& pc : ParticleCombination::cache)
            if (ParticleCombination::equivDown(PC, pc.lock()))
                SymmetrizationIndices_.push_back(pc.lock());
}

}

