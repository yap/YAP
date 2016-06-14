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
    Daughters_(daughters),
    DecayingParticle_(nullptr)
{
    // check daughter size
    if (Daughters_.empty())
        throw exceptions::Exception("No daughters", "DecayChannel::DecayChannel");
    if (Daughters_.size() == 1)
        throw exceptions::Exception("Only one daughter", "DecayChannel::DecayChannel");
    if (Daughters_.size() > 2)
        throw exceptions::Exception("More than two daughters", "DecayChannel::DecayChannel");

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
        // get daughter's list of particle combinations
        ParticleCombinationVector v_d = d->particleCombinations();

        // create vector copy them into
        ParticleCombinationVector v;
        v.reserve(v_d.size());
        
        // copy them if they aren't nullptr, aren't empty,
        // and aren't already contained in v, using equalDown to compare regardless of parent.
        std::copy_if(v_d.begin(), v_d.end(), std::back_inserter(v),
                     [&](const ParticleCombinationVector::value_type & pc)
                     {return pc and !pc->indices().empty()
                             and std::none_of(v.begin(), v.end(), std::bind(ParticleCombination::equalDown, std::placeholders::_1, pc));});

        if (v.empty())
            throw exceptions::Exception(std::string("No ParticleCombinations for daughter ") + to_string(*d)
                                        + " in DecayChannel " + to_string(*this),
                                        "DecayChannel::DecayChannel");

        PCs.push_back(v);
    }

    // create ParticleCombination's of parent
    /// \todo remove hardcoding for two daughters so applies to n daughters?
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
                if (!b_a.expired() and any_of(particleCombinations(), b_a.lock(), ParticleCombination::equalBySharedPointer))
                    // if (B,A) already added, don't proceed to adding (A,B)
                    continue;
            }

            // create (A,B), ParticleCombinationCache::composite copies PCA and PCB,
            // setting the parents of both to the newly created ParticleCombination
            addParticleCombination(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite({PCA, PCB}));
        }
    }
}

//-------------------------
void DecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // if pc already possessed, do nothing
    if (any_of(particleCombinations(), pc, ParticleCombination::equalBySharedPointer))
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

    // add to SpinAmplitude's
    for (auto& sa : SpinAmplitudes_)
        sa->addParticleCombination(pc);
}

//-------------------------
void DecayChannel::pruneParticleCombinations()
{
    if (!model())
        throw exceptions::Exception("Model not set", "DecayChannel::pruneParticleCombinations");

    // remove entries that don't trace back to the ISP
    for (auto it = ParticleCombinations_.begin(); it != ParticleCombinations_.end(); ) {
        // find the top-most parent
        auto pc = *it;
        while (pc->parent())
            pc = pc->parent();
        // check if it's not an ISP
        if (pc->indices().size() != model()->finalStateParticles().size())
            // erase
            it = ParticleCombinations_.erase(it);
        else
            it++;
    }

    if (ParticleCombinations_.empty())
        throw exceptions::Exception("ParticleCombinations empty after pruning", "DecayChannel::pruneParticleCombinations");

    for (auto& p : Daughters_)
        p->pruneParticleCombinations();
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

    // if SpinAmplitude's have already been added by hand, don't add automatically
    if (SpinAmplitudes_.empty()) {

        auto two_J = DecayingParticle_->quantumNumbers().twoJ();
        auto two_j1 = Daughters_[0]->quantumNumbers().twoJ();
        auto two_j2 = Daughters_[1]->quantumNumbers().twoJ();

        // create spin amplitudes
        // loop over possible S: |j1-j2| <= S <= (j1+j2)
        for (unsigned two_S = std::abs<int>(two_j1 - two_j2); two_S <= two_j1 + two_j2; two_S += 2) {
            // loop over possible L: |J-s| <= L <= (J+s)
            for (unsigned L = std::abs<int>(two_J - two_S) / 2; L <= (two_J + two_S) / 2; ++L) {
                // add SpinAmplitude retrieved from cache
                addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(two_J, two_j1, two_j2, L, two_S));
            }
        }
    } else {

        // check DecayingPartcle_'s quantum numbers against existing SpinAmplitude's
        if (DecayingParticle_->quantumNumbers().twoJ() != SpinAmplitudes_[0]->initialTwoJ())
            throw exceptions::Exception("Spins don't match ", "DecayChannel::setDecayingParticle");

    }

    // let DecayingParticle know to create a BlattWeisskopf objects for necessary orbital angular momenta
    for (auto& sa : SpinAmplitudes_)
        DecayingParticle_->storeBlattWeisskopf(sa->L());
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
void DecayChannel::addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa)
{
    // check number of daughters
    if (sa->finalTwoJ().size() != Daughters_.size())
        throw exceptions::Exception("Number of daughters doesn't match", "DecayChannel::addSpinAmplitude");

    /// \todo quantum numbers more completely?
    // check against daughter quantum numbers
    for (size_t i = 0; i < Daughters_.size(); ++i)
        if (Daughters_[i]->quantumNumbers().twoJ() != sa->finalTwoJ()[i])
            throw exceptions::Exception("Spins don't match daughter's", "DecayChannel::addSpinAmplitude");

    // check against DecayingParticle_ if set
    if (DecayingParticle_) {
        if (DecayingParticle_->quantumNumbers().twoJ() != sa->initialTwoJ())
            throw exceptions::Exception("Spins don't match DecayingParticle", "DecayChannel::addSpinAmplitude");
    } else {
        // else check against previously added SpinAmplitude's initial quantum numbers
        if (!SpinAmplitudes_.empty() and SpinAmplitudes_[0]->initialTwoJ() != sa->initialTwoJ())
            throw exceptions::Exception("Spins don't match previously added", "DecayChannel::addSpinAmplitude");
    }

    // add to SpinAmplitudes_
    SpinAmplitudes_.push_back(sa);

    // add this' ParticleCombination's to it
    for (auto& pc : particleCombinations())
        SpinAmplitudes_.back()->addParticleCombination(pc);
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

    // loop over SpinAmplitude's
    for (const auto& sa : SpinAmplitudes_) {
        // check SpinAmplitude
        if (!sa) {
            FLOG(ERROR) << "A SpinAmplitude is empty";
            C &= false;
        } else
            C &= sa->consistent();
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
    // auto& saV = dc.spinAmplitudes();
    // if (saV.empty())
    //     return s;
    // s += " " + to_string(saV);
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
