#include "DecayChannel.h"

#include "container_utils.h"
#include "CachedDataValue.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "Particle.h"
#include "ParticleCombinationCache.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <assert.h>
#include <stdexcept>

namespace yap {

//-------------------------
DecayChannel::AmplitudePair::AmplitudePair(DecayChannel* dc, std::complex<double> free) :
    Fixed(ComplexCachedDataValue::create(dc)),
    Free(std::make_shared<ComplexParameter>(free))
{
}

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

    // check that first daughter's ISP is not nullptr
    if (Daughters_[0]->initialStateParticle() == nullptr)
        throw exceptions::Exception(std::string("InitialStateParticle unset in ") + to_string(*Daughters_[0]),
                                    "DecayChannel::DecayChannel");
    // check that all daughters have same ISP (trivially checks first daughter against itself)
    for (auto& d : Daughters_)
        if (d->initialStateParticle() != Daughters_[0]->initialStateParticle())
            throw exceptions::Exception("InitialStateParticle mismatch", "DecayChannel::DecayChannel");

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
                auto b_a = initialStateParticle()->particleCombinationCache().find({PCB, PCA});
                // if b_a is not in cache, it can't be in SymmetrizationIndices_
                if (!b_a.expired() and hasParticleCombination(b_a.lock()))
                    // if (B,A) already added, don't proceed to adding (A,B)
                    continue;
            }

            // create (A,B), ParticleCombinationCache::composite copies PCA and PCB,
            // setting the parents of both to the newly created ParticleCombination
            addParticleCombination(initialStateParticle()->particleCombinationCache().composite({PCA, PCB}));
        }
    }

    // register with ISP
    addToInitialStateParticle();
}

//-------------------------
void DecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    DataAccessor::addParticleCombination(pc);

    // add to pc's daughters to daughter particles;
    // pc's daughters have their parents set correctly.
    for (size_t i = 0; i < pc->daughters().size(); ++i)
        Daughters_[i]->addParticleCombination(pc->daughters()[i]);

    // add to SpinAmplitude's (keys of Amplitudes_)
    for (auto& kv : Amplitudes_)
        kv.first->addParticleCombination(pc);
}

//-------------------------
void DecayChannel::setDecayingParticle(DecayingParticle* dp)
{
    DecayingParticle_ = dp;
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle is nullptr", "DecayChannel::setDecayingParticle");

    // if SpinAmplitude's have already been added by hand, don't add automatically
    if (!Amplitudes_.empty())
        return;

    auto iQ = DecayingParticle_->quantumNumbers();
    auto d1Q = Daughters_[0]->quantumNumbers();
    auto d2Q = Daughters_[1]->quantumNumbers();

    // create spin amplitudes
    // loop over possible S: |j1-j2| <= S <= (j1+j2)
    for (unsigned two_S = std::abs<int>(d1Q.twoJ() - d2Q.twoJ()); two_S <= d1Q.twoJ() + d2Q.twoJ(); two_S += 2)
        // loop over possible L: |J-s| <= L <= (J+s)
        for (unsigned L = std::abs<int>(iQ.twoJ() - two_S) / 2; L <= (iQ.twoJ() + two_S) / 2; ++L)
            // add SpinAmplitude retrieved from cache
            addSpinAmplitude(initialStateParticle()->spinAmplitudeCache().spinAmplitude(iQ, d1Q, d2Q, L, two_S));
}

//-------------------------
void DecayChannel::addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa)
{
    // check number of daughters
    if (sa->finalQuantumNumbers().size() != Daughters_.size())
        throw exceptions::Exception("Number of daughters doesn't match", "DecayChannel::addSpinAmplitude");

    // check against daughter quantum numbers
    for (size_t i = 0; i < Daughters_.size(); ++i)
        if (Daughters_[i]->quantumNumbers() != sa->finalQuantumNumbers()[i])
            throw exceptions::Exception(std::string("QuantumNumbers don't match daughter ") + std::to_string(i),
                                        "DecayChannel::addSpinAmplitude");

    // check against DecayingParticle_ if set
    if (DecayingParticle_) {
        if (DecayingParticle_->quantumNumbers() != sa->initialQuantumNumbers())
            throw exceptions::Exception("QuantumNumbers don't match DecayingParticle", "DecayChannel::addSpinAmplitude");
    } else {
        // else check against previously added SpinAmplitude's initial quantum numbers
        if (!Amplitudes_.empty() and Amplitudes_.begin()->first->initialQuantumNumbers() != sa->initialQuantumNumbers())
            throw exceptions::Exception("QuantumNumbers don't match previously added", "DecayChannel::addSpinAmplitude");
    }

    // add this' ParticleCombination's to it
    for (auto& pc : particleCombinations())
        sa -> addParticleCombination(pc);

    // create vector of amplitude pairs, one for each spin projection in the SpinAmplitude
    AmplitudePairMap apM;
    for (auto& two_m : sa->twoM()) {

        // add TotalAmplitude for two_m if needed
        if (TotalAmplitudes_.find(two_m) == TotalAmplitudes_.end())
            TotalAmplitudes_[two_m] = ComplexCachedDataValue::create(this);

        // insert new AmplitudePair into map, retaining the added Amplitude pair
        auto ap = (apM.insert(std::make_pair(two_m, std::move(AmplitudePair(this))))).first->second;

        /// add SpinAmplitude's cached amplitudes as a dependency for the Fixed amplitude
        ap.Fixed->addDependencies(sa->amplitudeSet());

        // add daughter dependencies to the fixed amplitude
        for (int i = 0; i < int(Daughters_.size()); ++i)
            if (auto d = std::dynamic_pointer_cast<DecayingParticle>(Daughters_[i]))
                for (auto& c : d->CachedDataValuesItDependsOn())
                    ap.Fixed->addDependency(c, i);

        // add to TotalAmplitudes_[two_m]'s dependencies
        TotalAmplitudes_[two_m]->addDependency(ap.Free);
        TotalAmplitudes_[two_m]->addDependency(ap.Fixed);
    }

    // add to Amplitudes_
    Amplitudes_.insert(std::make_pair(sa, std::move(apM)));
}

//-------------------------
ComplexParameterVector DecayChannel::freeAmplitudes()
{
    ComplexParameterVector V;
    for (auto& apM_kv : Amplitudes_)
        for (auto& ap_kv : apM_kv.second)
            V.push_back(ap_kv.second.Free);
    return V;
}

//-------------------------
std::complex<double> DecayChannel::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
        int two_m, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    auto& totAmp = TotalAmplitudes_.at(two_m);

    if (totAmp->calculationStatus(pc, symIndex, dataPartitionIndex) != kUncalculated) {
        FLOG(DEBUG) << "\nused cached fixed amplitude for " << *this << " for " << *pc << " : " << totAmp->value(d, symIndex);
        return totAmp->value(d, symIndex);
    }

    // sum over L-S combinations
    // LOOP_0 = sum_{L, S} BlattWeisskopf_L * free_amp(L, S, m) * LOOP_1
    // kvA = pair< shared_ptr<SpinAmplitude>, vector<AmplitudePair> >
    std::complex<double> A = Complex_0;
    for (auto& kvA : Amplitudes_) {

        auto ap = kvA.second.at(two_m); // AmplitudePair for spin projection m

        // if already calculated

        if (ap.Fixed->calculationStatus(pc, symIndex, dataPartitionIndex) != kUncalculated) {
            // Fixed is already calculated, simply retrieve from cache
            A += ap.Free->value() * ap.Fixed->value(d, symIndex);
            // DEBUG("DecayChannel::amplitude - use cached fixed amplitude for " << *this << " " << *pc << " = " << ap.Fixed->value(d, symIndex));
            continue;
        }

        // else calculate

        auto sa = kvA.first; // SpinAmplitude

        // get map of SpinProjectionPair's to cached spin amplitudes
        const auto& m = sa->amplitudes().at(two_m);
        auto sa_symIndex = sa->symmetrizationIndex(pc);

        // sum over daughter spin projection combinations (m1, m2)
        // LOOP_1 = sum_{m1, m2} SpinAmplitude_{L, S, m, m1, m2}(d) * amp_{daughter1}(m1) * amp_{daughter2}(m2)
        // kvM = pair <SpinProjectionPair, ComplexCachedDataValue>
        std::complex<double> a = Complex_0;
        for (auto& kvM : m) {
            // retrieve cached spin amplitude from data point
            auto amp = kvM.second->value(d, sa_symIndex);

            // loop over daughters, multiplying by their amplitudes for their spin projections
            const auto& spp = kvM.first; // SpinProjectionPair
            for (size_t i = 0; i < spp.size(); ++i)
                amp *= Daughters_[i]->amplitude(d, pc->daughters()[i], spp[i], dataPartitionIndex);

            // add into total amplitude thus far
            a += amp;
        }

        // multiply sum by Blatt-Weisskopf factor for orbital angular momentum L
        a *= DecayingParticle_->BlattWeisskopfs_[sa->L()]->amplitude(d, pc, dataPartitionIndex);

        // store result
        ap.Fixed->setValue(a, d, symIndex, dataPartitionIndex);

        // add into total amplitude A
        A += a * ap.Free->value();
    }

    // store result
    totAmp->setValue(A, d, symIndex, dataPartitionIndex);

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
            if (kv.first->finalQuantumNumbers().size() != Daughters_.size()) {
                FLOG(ERROR) << "quantum numbers object and daughters object size mismatch";
                C &= false;
            }

            // check if QuantumNumbers of SpinAmplitude objects match with Particles
            if (kv.first->initialQuantumNumbers() != decayingParticle()->quantumNumbers()) {
                FLOG(ERROR) << "quantum numbers of parent " << decayingParticle()->quantumNumbers()
                            << " and SpinAmplitude " << kv.first->initialQuantumNumbers() << " don't match.";
                C &= false;
            }

            for (size_t i = 0; i < Daughters_.size(); ++i) {
                if (kv.first->finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
                    FLOG(ERROR) << "quantum numbers of daughter " << i << " " << Daughters_[i]->quantumNumbers()
                                << " and SpinAmplitude " << kv.first->finalQuantumNumbers()[i] << " don't match.";
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
    double mass_sum = std::accumulate(Daughters_.begin(), Daughters_.end(), 0.,
    [](double m, std::shared_ptr<Particle> d) {return m + (d ? d->mass()->value() : 0);});
    if (mass_sum > decayingParticle()->mass()->value()) {
        FLOG(ERROR) << "sum of daughter's masses (" << mass_sum << ")"
                    << "is bigger than resonance mass (" << decayingParticle()->mass()->value() << ").";
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
CachedDataValueSet DecayChannel::CachedDataValuesItDependsOn()
{
    CachedDataValueSet S;
    for (auto& kv : TotalAmplitudes_)
        S.insert(kv.second);
    return S;
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
    auto& saV = dc.spinAmplitudes();
    if (saV.empty())
        return s;
    s += " " + to_string(saV);
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

//-------------------------
/*CachedDataValuePcIndexSet DecayChannel::CachedDataValuesItDependsOn()
{
    CachedDataValuePcIndexSet set;
    set.insert(std::make_pair(FixedAmplitude_, -1));

    for (int i=0; i<int(Daughters_.size()); ++i) {
        auto daugh = std::dynamic_pointer_cast<DecayingParticle>(Daughters_[i]);
        if (!daugh)
            continue;
        for (auto& c : daugh->CachedDataValuesItDependsOn()) {
            if (c.second >= 0)
                LOG(FATAL) << "fatal error";
            set.insert(std::make_pair(c.first, i));
        }
    }

    return set;
}*/


}
