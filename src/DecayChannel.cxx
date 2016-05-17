#include "DecayChannel.h"

#include "BlattWeisskopf.h"
#include "container_utils.h"
#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
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
    DataAccessor(),
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
                if (!b_a.expired() and hasParticleCombination(b_a.lock()))
                    // if (B,A) already added, don't proceed to adding (A,B)
                    continue;
            }

            // create (A,B), ParticleCombinationCache::composite copies PCA and PCB,
            // setting the parents of both to the newly created ParticleCombination
            addParticleCombination(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite({PCA, PCB}));
        }
    }

    // register with Model
    addToModel();
}

//-------------------------
void DecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    DataAccessor::addParticleCombination(pc);

    // add pc's daughters to daughter particles;
    // pc's daughters have their parents set correctly.
    for (size_t i = 0; i < pc->daughters().size(); ++i)
        Daughters_[i]->addParticleCombination(pc->daughters()[i]);

    // add to SpinAmplitude's (keys of Amplitudes_)
    for (auto& kv : Amplitudes_)
        kv.first->addParticleCombination(pc);
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
    if (Amplitudes_.empty()) {

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
    }

    // let DecayingParticle know to create a BlattWeisskopf objects for necessary orbital angular momenta
    // loop over mapping of (spin amplitude) -> (amplitude pair map)
    for (auto& sa_apm : Amplitudes_)
        DecayingParticle_->storeBlattWeisskopf(sa_apm.first->L());
}

//-------------------------
const Model* DecayChannel::model() const
{
    return Daughters_[0]->model();
}

//-------------------------
void DecayChannel::addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa)
{
    // check number of daughters
    if (sa->finalTwoJ().size() != Daughters_.size())
        throw exceptions::Exception("Number of daughters doesn't match", "DecayChannel::addSpinAmplitude");

    /// \todo see what needs to be checked
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
        if (!Amplitudes_.empty() and Amplitudes_.begin()->first->initialTwoJ() != sa->initialTwoJ())
            throw exceptions::Exception("Spins don't match previously added", "DecayChannel::addSpinAmplitude");
    }

    // add this' ParticleCombination's to it
    for (auto& pc : particleCombinations())
        sa -> addParticleCombination(pc);

    // create new map_type::mapped_type (map of spin projection -> free amplitude)
    map_type::mapped_type m_fa;
    // create free amplitude for each spin projection
    for (auto& two_m : sa->twoM())
        m_fa.emplace(two_m, std::make_shared<ComplexParameter>(Complex_1));

    // add to Amplitudes_
    Amplitudes_.emplace(sa, std::move(m_fa));

}

//-------------------------
std::shared_ptr<ComplexParameter> DecayChannel::freeAmplitude(int two_M, unsigned l, unsigned two_s)
{
    auto it1 = std::find_if(Amplitudes_.begin(), Amplitudes_.end(), [&](const map_type::value_type & kv) {return kv.first->L() == l and kv.first->twoS() == two_s;});
    if (it1 == Amplitudes_.end())
        throw exceptions::Exception("No amplitude found for L = " + std::to_string(l) + " and S = " + spin_to_string(two_s), "DecayChannel::freeAmplitude");

    auto it2 = it1->second.find(two_M);
    if (it2 == it1->second.end())
        throw exceptions::Exception("No amplitude found for spin projection = " + spin_to_string(two_M), "DecayChannel::freeAmplitude");

    return it2->second;
}

//-------------------------
ComplexParameterVector DecayChannel::freeAmplitudes()
{
    ComplexParameterVector V;
    for (auto& sa_mfa : Amplitudes_)
        for (auto& m_fa : sa_mfa.second)
            V.push_back(m_fa.second);
    return V;
}

//-------------------------
std::complex<double> DecayChannel::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const
{
    std::complex<double> A = Complex_0;

    // sum over L-S combinations
    // LOOP_0 = sum_{L, S} BlattWeisskopf_L * free_amp(L, S, m) * LOOP_1
    // sa_mfa = pair< shared_ptr<SpinAmplitude>, map<spin projection, free amplitude> >
    for (const auto& sa_mfa : Amplitudes_) {

        const auto& sa = sa_mfa.first; // SpinAmplitude
        const auto& fa = sa_mfa.second.at(two_m); // free amplitude for spin projection m

        // get map of SpinProjectionPair's to cached spin amplitudes
        const auto& m_cdv_map = sa->amplitudes().at(two_m);
        const auto sa_symIndex = sa->symmetrizationIndex(pc);

        // sum over daughter spin projection combinations (m1, m2)
        // LOOP_1 = sum_{m1, m2} SpinAmplitude_{L, S, m, m1, m2}(d) * amp_{daughter1}(m1) * amp_{daughter2}(m2)
        // m_cdv = pair <SpinProjectionPair, ComplexCachedDataValue>
        std::complex<double> a = Complex_0;
        for (const auto& m_cdv : m_cdv_map) {
            // retrieve cached spin amplitude from data point
            auto amp = m_cdv.second->value(d, sa_symIndex);

            // loop over daughters, multiplying by their amplitudes for their spin projections
            const auto& spp = m_cdv.first; // SpinProjectionPair
            for (size_t i = 0; i < spp.size(); ++i)
                amp *= Daughters_[i]->amplitude(d, pc->daughters()[i], spp[i], sm);

            // add into total amplitude thus far
            a += amp;
        }

        // multiply sum by Blatt-Weisskopf factor for orbital angular momentum L
        a *= DecayingParticle_->BlattWeisskopfs_[sa->L()]->amplitude(d, pc, sm);

        // add into total amplitude A
        A += a * fa->value();

    }

    // and return it
    return A;
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool C = DataAccessor::consistent();

    // check number of daughters
    if (Daughters_.size() < 2) {
        FLOG(ERROR) << "invalid number of daughters (" << Daughters_.size() << " < 2).";
        C &= false;
    }

    // currently only allowing exactly two daughters
    if (Daughters_.size() > 2) {
        FLOG(ERROR) << "invalid number of daughters (" << Daughters_.size() << " > 2).";
        C &= false;
    }

    // compare number of daughters with sizes of ParticleCombinations
    auto pcs = particleCombinations();
    if (std::any_of(pcs.begin(), pcs.end(),
    [&](ParticleCombinationVector::value_type pc) {return pc->daughters().size() != Daughters_.size();})) {
        FLOG(ERROR) << "DecayChannel and its particleCombinations do not have the same number of daughters.";
        C &= false;
    }

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

    // loop over SpinAmplitude's, which are keys (first) to Amplitudes_ map
    for (auto& kv : Amplitudes_) {
        // check SpinAmplitude
        if (!kv.first) {
            FLOG(ERROR) << "A SpinAmplitude is empty";
            C &= false;
        } else {
            C &= kv.first->consistent();

            // check size of SpinAmplitude's quantum numbers against size of daughters
            if (kv.first->finalTwoJ().size() != Daughters_.size()) {
                FLOG(ERROR) << "quantum numbers object and daughters object size mismatch";
                C &= false;
            }

            // check if QuantumNumbers of SpinAmplitude objects match with Particles
            if (kv.first->initialTwoJ() != decayingParticle()->quantumNumbers().twoJ()) {
                FLOG(ERROR) << "spins of parent " << decayingParticle()->quantumNumbers().twoJ()
                            << " and SpinAmplitude " << kv.first->initialTwoJ() << " don't match.";
                C &= false;
            }

            for (size_t i = 0; i < Daughters_.size(); ++i) {
                if (kv.first->finalTwoJ()[i] != Daughters_[i]->quantumNumbers().twoJ()) {
                    FLOG(ERROR) << "spins of daughter " << i << " " << Daughters_[i]->quantumNumbers().twoJ()
                                << " and SpinAmplitude " << kv.first->finalTwoJ()[i] << " don't match.";
                    C &= false;
                }
            }

            // check its corresponding BlattWeisskopf
            if (DecayingParticle_) {
                // look for BlattWeisskopf object with corresponding orbital angular momentum
                auto it = DecayingParticle_->BlattWeisskopfs_.find(kv.first->L());
                if (it == DecayingParticle_->BlattWeisskopfs_.end() or !it->second) {
                    FLOG(ERROR) << "Could not find BlattWeisskopf object with L = " << kv.first->L();
                    C &= false;
                } else {
                    C &= it->second->consistent();
                }
            }
        }
    }

    // check masses
    // double mass_sum = std::accumulate(Daughters_.begin(), Daughters_.end(), 0.,
    // [](double m, std::shared_ptr<Particle> d) {return m + (d ? d->mass()->value() : 0);});
    // if (mass_sum > decayingParticle()->mass()->value()) {
    //     FLOG(ERROR) << "sum of daughter's masses (" << mass_sum << ")"
    //                 << "is bigger than decaying particle's mass (" << decayingParticle()->mass()->value() << ").";
    //     C &= false;
    // }

    // check charge conservation
    int daughtersQ(0);
    for (auto& d : Daughters_)
        daughtersQ += d->quantumNumbers().Q();
    if (DecayingParticle_->quantumNumbers().Q() != daughtersQ) {
        FLOG(ERROR) << "charge conservation violated";
        C &= false;
    }

    return C;
}

//-------------------------
SpinAmplitudeVector DecayChannel::spinAmplitudes()
{
    SpinAmplitudeVector V;
    for (auto& kv : Amplitudes_)
        V.push_back(kv.first);
    return V;
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
