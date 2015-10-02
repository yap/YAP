#include "FinalStateParticle.h"

#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double mass, std::string name, std::vector<ParticleIndex>& indices)
    : Particle(q, mass, name)
{
    for (ParticleIndex i : indices)
        addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(i));
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
void FinalStateParticle::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    // check if not yet there
    if (std::find(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), c) == SymmetrizationIndices_.end()) {
        SymmetrizationIndices_.push_back(c);
    } else {
        LOG(WARNING) << "FinalStateParticle::addSymmetrizationIndex() - Index already existing!";
    }
}

//-------------------------
void FinalStateParticle::setSymmetrizationIndexParents()
{
    std::vector<std::shared_ptr<const ParticleCombination> > PCs = SymmetrizationIndices_;

    // check if already set
    if (PCs[0]->parent() != nullptr)
        return;

    SymmetrizationIndices_.clear();

    for (auto& PC : PCs) {
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivDown(PC, pc)) {
                //std::cout << "  add " << std::string(*pc) << " to fsp " << name() << "\n";
                SymmetrizationIndices_.push_back(pc);
            }
        }
    }

}

}

