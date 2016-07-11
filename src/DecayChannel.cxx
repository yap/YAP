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
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"
#include "Spin.h"
#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"
#include "StatusManager.h"

#include <algorithm>
#include <iterator>

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

    // check no Daughters_ are empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), std::logical_not<ParticleVector::value_type>()))
        throw exceptions::Exception("Empty daughter", "DecayChannel::DecayChannel");

    // if more than two daughters, check all are spin 0
    if (Daughters_.size() > 2 and
        std::any_of(Daughters_.begin(), Daughters_.end(), [](const ParticleVector::value_type & d) {return d->quantumNumbers().twoJ() != 0;})) {
        LOG(ERROR) << "to create so-called \"nonresonant\" decays with nonzero-spin daughters, "
                   << "create a DecayingParticle for subcontent and use two-particle decays.";
        throw exceptions::Exception("Attempted to create nonresonant decay with spinfull daughters", "DecayChannel::DecayChannel");
    }

    // check that (first daughter's) Model is not nullptr
    if (model() == nullptr)
        throw exceptions::Exception(std::string("Model unset in ") + to_string(*Daughters_[0]),
                                    "DecayChannel::DecayChannel");

    // check that all daughters have same Model (trivially checks first daughter against itself)
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [&](const ParticleVector::value_type & d) {return model() != d->model();}))
        throw exceptions::Exception("Model mismatch", "DecayChannel::DecayChannel");

    //////////////////////////////////////////////////
    // create ParticleCombination's for parent
    //

    // create vector of unique entries of daughters
    auto u_daughters = Daughters_;
    std::sort(u_daughters.begin(), u_daughters.end());
    u_daughters.erase(std::unique(u_daughters.begin(), u_daughters.end()), u_daughters.end());

    // collect combinations of ParticleCombination's of daughters
    // map: daughter -> vector of combinations of particle combinations
    std::map<ParticleVector::value_type, std::vector<ParticleCombinationVector> > d_CV_map;
    std::vector<ParticleCombinationVector> PCs;
    for (auto d : u_daughters) {
        // get daughter's list of particle combinations
        ParticleCombinationVector v_d = d->particleCombinations();

        // create vector to copy them into
        ParticleCombinationVector v;
        v.reserve(v_d.size());

        // copy them if they aren't nullptr, aren't empty,
        // and aren't already contained in v, using equalDown to compare regardless of parent.
        std::copy_if(v_d.begin(), v_d.end(), std::back_inserter(v),
                     [&](const ParticleCombinationVector::value_type & pc) {
                         return pc and !pc->indices().empty() and
                             std::none_of(v.begin(), v.end(), std::bind(equal_down, std::placeholders::_1, pc));
                     });

        if (v.empty())
            throw exceptions::Exception("No ParticleCombinations for daughter " + to_string(*d)
                                        + " in DecayChannel " + to_string(*this),
                                        "DecayChannel::DecayChannel");

        // create combinations of ParticleCombinations to accomodate occurance of daughter in final state
        auto cv = combinations(v, std::count(Daughters_.begin(), Daughters_.end(), d));
        if (cv.empty())
            throw exceptions::Exception("Could not form combinations of ParticleCombinations for daughter " + to_string(*d),
                                        "DecayChannel::DecayChannel");

        // place combinations into map:
        d_CV_map.emplace(d, cv);
    }

    // create vector of iterators to combinations for each unique daughter
    // initialized to begin()
    std::vector<std::vector<ParticleCombinationVector>::const_iterator> Its;
    Its.reserve(u_daughters.size());
    std::transform(u_daughters.begin(), u_daughters.end(), std::back_inserter(Its),
    [&](const ParticleVector::value_type & ud) {return d_CV_map[ud].begin();});

    // use 'odometer" style increment of iterators:
    // repeat until last iterator hits its end()
    while (Its.back() != d_CV_map[u_daughters.back()].end()) {

        // fill vector of ParticleCombination's to use to create parent ParticleCombination from
        ParticleCombinationVector pcv(Daughters_.size());
        // loop through iterators to combinations
        for  (size_t d = 0; d < Its.size(); ++d) {
            auto prev_it = Daughters_.begin();
            // loop over elements in particular combination
            for (const auto& pc : *Its[d]) {
                // find daughter to set value for, starting from one beyond previously found daughter
                auto it = std::find(prev_it, Daughters_.end(), u_daughters[d]);
                if (it == Daughters_.end())
                    throw exceptions::Exception("failed to create particle combinations", "DecayChannel::DecayChannel");
                // set daughter's ParticleCombination
                pcv[it - Daughters_.begin()] = pc;
                // increment prev_it to one beyond last found daughter
                prev_it = it + 1;
            }
        }

        // check vector for possible overlaps
        if (disjoint(pcv))
            // create new ParticleCombination for parent;
            // ParticleCombinationCache::composite copies all elts of
            // pcv setting the parents of all copies to the newly
            // created ParticleCombination
            addParticleCombination(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite(pcv));

        // increment the "odometer":
        // increment first iterator
        ++Its[0];
        // loop through iterators: if current is at end(), reset it and increment next one
        for (size_t i = 0; (i < Its.size() - 1) and Its[i] == d_CV_map[u_daughters[i]].end(); ++i) {
            Its[i] = d_CV_map[u_daughters[i]].begin();
            ++Its[i + 1];
        }
        // the "odometer" has ticked over completely if the last iterator is now at its end()
    }

    // check that at least one combination was created above
    if (particleCombinations().empty())
        throw exceptions::Exception("ParticleCombinations_ is empty", "DecayChannel::DecayChannel");

    // if more than two daughters, add nonresonant spin amplitude
    if (Daughters_.size() > 2)
        // add SpinAmplitude retrieved from cache
        addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(0, spins(Daughters_), 0, 0));
}

//-------------------------
void DecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // if pc already possessed, do nothing
    if (any_of(particleCombinations(), pc))
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
    prune_particle_combinations(ParticleCombinations_);

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
    if (DecayingParticle_)
        throw exceptions::Exception("DecayingParticle is already set", "DecayChannel::setDecayingParticle");

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
        auto two_j = spins(Daughters_);

        // create spin amplitudes
        // loop over possible S: |j1-j2| <= S <= (j1+j2)
        for (unsigned two_S = std::abs<int>(two_j[0] - two_j[1]); two_S <= two_j[0] + two_j[1]; two_S += 2) {
            // loop over possible L: |J-s| <= L <= (J+s)
            for (unsigned L = std::abs<int>(two_J - two_S) / 2; L <= (two_J + two_S) / 2; ++L) {
                // add SpinAmplitude retrieved from cache
                addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(two_J, two_j, L, two_S));
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
    // if spin amplitude already held, do nothing
    if (std::find(SpinAmplitudes_.begin(), SpinAmplitudes_.end(), sa) != SpinAmplitudes_.end())
        return;

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
FinalStateParticleVector DecayChannel::finalStateParticles() const
{
    FinalStateParticleVector fsps;

    for (std::shared_ptr<Particle> d : Daughters_) {

        // if daughter is fsp
        if (std::dynamic_pointer_cast<FinalStateParticle>(d)) {
            fsps.push_back(std::static_pointer_cast<FinalStateParticle>(d));

        } else if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
            auto ddaughters = std::dynamic_pointer_cast<DecayingParticle>(d)->finalStateParticles(0);
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());

        } else {
            throw exceptions::Exception("Invalid daughter", "DecayChannel::DecayChannel");
        }

    }

    return fsps;
}

}
