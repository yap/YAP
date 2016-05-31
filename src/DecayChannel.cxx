#include "DecayChannel.h"

#include "BlattWeisskopf.h"
#include "container_utils.h"
#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombinationCache.h"
#include "spin.h"
#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(const ParticleVector& daughters) :
    ReportsParticleCombinations(),
    Daughters_(daughters),
    DecayingParticle_(nullptr)
{
    // check daughter size
    if (Daughters_.empty())
        throw exceptions::Exception("No daughters", "DecayChannel::DecayChannel");
    if (Daughters_.size() == 1)
        throw exceptions::Exception("Only one daughter", "DecayChannel::DecayChannel");

    // check no Daughters_ are empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [](std::shared_ptr<Particle> d) {return !d;}))
    throw exceptions::Exception("Empty daughter", "DecayChannel::DecayChannel");

    // check that first daughter's Model is not nullptr
    if (Daughters_[0]->model() == nullptr)
        throw exceptions::Exception(std::string("Model unset in ") + to_string(*Daughters_[0]),
                                    "DecayChannel::DecayChannel");

    // check that all daughters have same Model (trivially checks first daughter against itself)
    for (auto& d : Daughters_)
        if (d->model() != Daughters_[0]->model())
            throw exceptions::Exception("Model mismatch", "DecayChannel::DecayChannel");

    // collect ParticleCombination's of daughters
    std::vector<ParticleCombinationVector> PCs;
    for (auto d : Daughters_) {
        ParticleCombinationVector v;
        ParticleCombinationVector v_d = d->particleCombinations();
        for (auto pc : v_d) {
            // check for empty indices
            if (pc->indices().empty())
                throw exceptions::Exception("ParticleCombination has empty indices", "DecayChannel::DecayChannel");
            // ignore PC's that differ only by parent from ones already accounted for
            if (std::none_of(v.begin(), v.end(), [&](const std::shared_ptr<ParticleCombination>& A) {return ParticleCombination::equivDown(A, pc);}))
            v.push_back(pc);
        }
        if (v.empty())
            throw exceptions::Exception(std::string("No ParticleCombinations for daughter ") + to_string(*d)
                                        + " in DecayChannel " + to_string(*this),
                                        "DecayChannel::DecayChannel");
        PCs.push_back(v);
    }

    // create ParticleCombination's of parent
    if (PCs.size() == 2) {
        for (auto& PCA : PCs[0]) {
            for (auto& PCB : PCs[1]) {

                // check that PCA and PCB don't overlap in FSP content
                if (overlap(PCA->indices(), PCB->indices()))
                    continue;

                // for identical particles, check if swapped particle combination is already added
                if (Daughters_[0] == Daughters_[1]) {
                    // get (B,A) combination from cache
                    auto b_a = model()->particleCombinationCache().find({PCB, PCA});
                    // if b_a is not in cache, it can't be in SymmetrizationIndices_
                    if (!b_a.expired() and hasParticleCombination(b_a.lock(), ParticleCombination::equivBySharedPointer))
                        // if (B,A) already added, don't proceed to adding (A,B)
                        continue;
                }

                // create (A,B), ParticleCombinationCache::composite copies PCA and PCB,
                // setting the parents of both to the newly created ParticleCombination
                addParticleCombination(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite({PCA, PCB}));
            }
        }
    }
    else { // > 2 daughters

        std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > PCsVec;
        for (size_t i = 0; i < PCs.size(); ++i) {
            std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > vec;
            for (auto& p : PCs[i])
                vec.push_back(std::make_pair(Daughters_[i], ParticleCombinationVector(1, p)));

            PCsVec.push_back(vec);
        }

        auto possibleParents = combinations(PCsVec, const_cast<Model*>(static_cast<const DecayChannel*>(this)->model()));

        for (auto& parents : possibleParents) {
            for (auto& p : parents) {
                DEBUG("parent: " << to_string(p.second));
                addParticleCombination(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite(p.second));
            }
        }
    } // end else // > 2 daughters
}

//-------------------------
void DecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // if pc already possessed, do nothing
    if (hasParticleCombination(pc, ParticleCombination::equivBySharedPointer))
        return;

    // check number of daughters in pc
    if (pc->daughters().size() != Daughters_.size())
        throw exceptions::Exception("ParticleCombination has wrong number of daughters ("
                                    + std::to_string(pc->daughters().size()) + " != "
                                    + std::to_string(Daughters_.size()) + ")",
                                    "DecayChannel::addParticleCombination");

    ParticleCombinations_.push_back(pc);

    // add pc's daughters to daughter particles;
    // pc's daughters have their parents set correctly.
    for (size_t i = 0; i < pc->daughters().size(); ++i)
        Daughters_[i]->addParticleCombination(pc->daughters()[i]);

}

//-------------------------
void DecayChannel::fixSolitaryFreeAmplitudes()
{
    for (auto& d : Daughters_)
        if (std::dynamic_pointer_cast<DecayingParticle>(d))
            std::static_pointer_cast<DecayingParticle>(d)->fixSolitaryFreeAmplitudes();
}

//-------------------------
void DecayChannel::setDecayingParticle(DecayingParticle* dp)
{
    DecayingParticle_ = dp;
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle is nullptr", "DecayChannel::setDecayingParticle");

    // check charge conservation
    int q = 0;
    for (const auto& d : Daughters_)
        q += d->quantumNumbers().Q();
    if (DecayingParticle_->quantumNumbers().Q() != q)
        throw exceptions::Exception("Charge not conserved: " + std::to_string(DecayingParticle_->quantumNumbers().Q()) + " != " + std::to_string(q)
                                    + " in " + to_string(*this),
                                    "DecayChannel::setDecayingParticle");

}

//-------------------------
const Model* DecayChannel::model() const
{
    return Daughters_[0]->model();
}

//-------------------------
FreeAmplitudeSet DecayChannel::freeAmplitudes() const
{
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle is nullptr", "DecayChannel::freeAmplitudes");

    // get set from decaying particle
    return find(DecayingParticle_->freeAmplitudes(), this);
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool C = true;

    // check no daughters is empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [](std::shared_ptr<Particle> d) {return !d;})) {
        FLOG(ERROR) << "null pointer in daughters vector.";
        C &= false;
    }
    // check daughters
    std::for_each(Daughters_.begin(), Daughters_.end(), [&](std::shared_ptr<Particle> d) {if (d) C &= d->consistent();});

    // check DecayingParticle_ is set
    if (!DecayingParticle_) {
        FLOG(ERROR) << "DecayingParticle is unset.";
        C &= false;
    }

    return C;
}

//-------------------------
std::string to_string(const DecayChannel& dc)
{
    std::string s = "(";
    if (dc.decayingParticle())
        s += dc.decayingParticle()->name() + " -> ";
    if (dc.daughters().empty())
        s += "[nothing]";
    else {
        for (auto& d : dc.daughters())
            s += d->name() + " ";
        s.erase(s.size() - 1, 1);
    }
    s += ")";

    return s;
}

//-------------------------
std::vector<std::shared_ptr<FinalStateParticle> > DecayChannel::finalStateParticles() const
{
    std::vector<std::shared_ptr<FinalStateParticle> > fsps;

    for (std::shared_ptr<Particle> d : Daughters_) {

        // if daughter is fsp
        if (std::dynamic_pointer_cast<FinalStateParticle>(d)) {
            fsps.push_back(std::static_pointer_cast<FinalStateParticle>(d));

        } else if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
            auto ddaughters = std::dynamic_pointer_cast<DecayingParticle>(d)->finalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());

        } else {
            FLOG(ERROR) << "Daughter is neither a FinalStateParticle nor a DecayingParticle. DecayChannel is inconsistent.";
            throw exceptions::Exception("Invalid daughter", "DecayChannel::DecayChannel");
        }
    }

    return fsps;
}

}
