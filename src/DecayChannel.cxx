#include "DecayChannel.h"

#include "container_utils.h"
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
DecayChannel::AmplitudePair::AmplitudePair(DecayChannel* dc, std::complex<double> fixed = Complex_1) :
    Fixed(std::shared_ptr<ComplexCachedDataValue>(dc)),
    Free(std::shared_ptr<ComplexParameter(free))
{
}

//-------------------------
DecayChannel::DecayChannel(ParticleVector daughters)
    DataAccessor(),
    Daughters_(daughters),
    BlattWeisskopf_(std::make_shared<BlattWeisskopf>(this)),
    SpinAmplitude_(spinAmplitude),
    FreeAmplitude_(std::make_shared<ComplexParameter>(Complex_1)),
    FixedAmplitude_(std::make_shared<ComplexCachedDataValue>(this)),
    DecayingParticle_(nullptr)
{
    // check SpinAmplitude
    if (!SpinAmplitude_)
        throw exceptions::Exception("SpinAmplitude unset", "DecayChannel::DecayChannel");
    // check that ISP is set (it's taken over from SpinAmplitude_)
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "DecayChannel::DecayChannel");

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
            // ignore PC's that differ only by parent or lambda from ones already accounted for
            if (std::none_of(v.begin(), v.end(), [&](const std::shared_ptr<ParticleCombination>& A) {return ParticleCombination::equivDownButLambda(A, pc);}))
            v.push_back(pc);
        }
        if (v.empty())
            throw exceptions::Exception(std::string("No ParticleCombinations for daughter ") + to_string(*d)
                                        + " in DecayChannel " + to_string(*this),
                                        "DecayChannel::DecayChannel");
        PCs.push_back(v);
    }

    // create ParticleCombnation's of parent
    /// \todo remove hardcoding for two daughters so applies to n daughters
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
                if (!b_a.expired() and hasSymmetrizationIndex(b_a.lock()))
                    // if (B,A) already added, don't proceed to adding (A,B)
                    continue;
            }

            // create (A,B), ParticleCombinationCache::composite copies PCA and PCB,
            // setting the parents of both to the newly created ParticleCombination
            auto a_b = initialStateParticle()->particleCombinationCache().composite({PCA, PCB});
            // add it to the spin amplitude, which returns a vector of ParticleCombinations with helicities set
            ParticleCombinationVector V = SpinAmplitude_->addSymmetrizationIndices(a_b);
            // add each resulting ParticleCombination to this
            for (auto v : V)
                addSymmetrizationIndex(v);
        }
    }
}

//-------------------------
void DecayChannel::setDecayingParticle(DecayingParticle* dp)
{
    DecayingParticle_ = dp;
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle is nullptr", "DecayChannel::setDecayingParticle");

    auto isp = Daughters_[0]->initialStateParticle();
    auto iQ = DecayingParticle_->quantumNumbers();
    auto d1Q = Daughters_[0]->quantumNumbers();
    auto d2Q = Daughters_[1]->quantumNumbers();
    
    // create spin amplitudes
    // loop over possible S: |j1-j2| <= S <= (j1+j2)
    for (unsigned two_S = std::abs<int>(d1Q.twoJ() - d2Q.twoJ()); two_S <= d1Q.twoJ() + d2.twoJ(); two_S += 2) {

        // loop over possible L: |J-s| <= L <= (J+s)
        for (unsigned L = std::abs<int>(iQ.twoJ() - two_S) / 2; L <= (iQ.twoJ() + two_S) / 2; ++L) {

            // retrieve SpinAmplitude from cache
            auto sa = isp->spinAmplitudeCache().spinAmplitude(iQ, d1Q, d2Q, L, two_S);

            // add to Amplitudes, recording kv pair in A
            auto A = *(Amplitudes_.insert(std::make_pair(sa, AmplitudePair(this))).first);

            /// add spin amplitude's Amplitude_ as a dependency for the Fixed amplitude
            A.second.Fixed->addDependency(A.first->amplitude());
            
            // add daughter dependencies to the fixed amplitude
            for (int i = 0; i < int(Daughters_.size()); ++i)
                if (auto d = std::dynamic_pointer_cast<DecayingParticle>(Daughters_[i]))
                    for (auto& c : d->CachedDataValuesItDependsOn())
                        A.second.Fixed->addDependency(c, i);
            
        }
    }
         
    
    
    
    BlattWeisskopf_->setDependencies();
}

    
//-------------------------
std::complex<double> DecayChannel::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    DEBUG("DecayChannel::amplitude - " << *this << " " << *pc);

    unsigned symIndex = symmetrizationIndex(pc);

    if (FixedAmplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        std::complex<double> a = BlattWeisskopf_->amplitude(d, pc, dataPartitionIndex) * SpinAmplitude_->amplitude(d, pc);

        auto& pcDaughters = pc->daughters();
        for (unsigned i = 0; i < Daughters_.size(); ++i)
            a *= Daughters_[i]->amplitude(d, pcDaughters.at(i), dataPartitionIndex);

        FixedAmplitude_->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("DecayChannel::amplitude - calculated fixed amplitude for " << this << " " << *pc << " = " << a);
        return FreeAmplitude_->value() * a;
    }

    DEBUG("DecayChannel::amplitude - use cached fixed amplitude for " << *this << " " << *pc << " = " << FixedAmplitude_->value(d, symIndex));
    return FreeAmplitude_->value() * FixedAmplitude_->value(d, symIndex);
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

    // Check Blatt-Weisskopf object
    if (!BlattWeisskopf_) {
        FLOG(ERROR) << "BlattWeisskopf is nullptr";
        C &= false;
    } else {
        C &= BlattWeisskopf_->consistent();

        // check if BlattWeisskopf points back to this DecayChannel
        if (BlattWeisskopf_->decayChannel() != this) {
            FLOG(ERROR) << "BlattWeisskopf does not point back to this DecayChannel.";
            C &=  false;
        }
    }

    // Check SpinAmplitude object
    if (!SpinAmplitude_) {
        FLOG(ERROR) << "no SpinAmplitude object set.";
        C &= false;
    } else {
        C &= SpinAmplitude_->consistent();

        // check size of spin amplitude quantum numbers and size of daughters
        if (SpinAmplitude_->finalQuantumNumbers().size() != Daughters_.size()) {
            FLOG(ERROR) << "quantum numbers object and daughters object size mismatch";
            C &= false;
        }

        // check if QuantumNumbers of SpinAmplitude objects match with Particles
        if (SpinAmplitude_->initialQuantumNumbers() != decayingParticle()->quantumNumbers()) {
            FLOG(ERROR) << "quantum numbers of parent " << decayingParticle()->quantumNumbers()
                        << " and SpinAmplitude " << SpinAmplitude_->initialQuantumNumbers() << " don't match.";
            C &= false;
        }

        for (size_t i = 0; i < Daughters_.size(); ++i) {
            if (SpinAmplitude_->finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
                FLOG(ERROR) << "quantum numbers of daughter " << i << " " << Daughters_[i]->quantumNumbers()
                            << " and SpinAmplitude " << SpinAmplitude_->finalQuantumNumbers()[i] << " don't match.";
                C &= false;
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
DecayChannel::AmplitudePair& DecayChannel::amplitudes(const std::shared_ptr<SpinAmplitude>& sa)
{
    if (Amplitudes_.find(sa) != Amplitudes_.end())
        return Amplitudes_[sa].second;
    return AmplitudePair()
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
    if (dc.spinAmplitude())
        s += " " + to_string(*dc.spinAmplitude());
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
void DecayChannel::clearSymmetrizationIndices()
{
    DataAccessor::clearSymmetrizationIndices();
    BlattWeisskopf_->clearSymmetrizationIndices();
    SpinAmplitude_->clearSymmetrizationIndices();
}

//-------------------------
void DecayChannel::setSymmetrizationIndexParents()
{
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "DecayChannel::DecayChannel");

    ParticleCombinationVector chPCs = particleCombinations();

    // clean up PCs without parents
    ParticleCombinationVector chPCsParents = particleCombinations();
    auto it = chPCsParents.begin();
    while (it != chPCsParents.end()) {
        if (not (*it)->parent()) {
            it = chPCsParents.erase(it);
        } else
            ++it;
    }
    clearSymmetrizationIndices();

    for (auto& pc : chPCsParents)
        addSymmetrizationIndex(pc);


    for (auto& chPC : chPCs) {
        for (auto& wpc : initialStateParticle()->particleCombinationCache()) {

            if (wpc.expired())
                continue;

            auto pc = wpc.lock();

            if (!ParticleCombination::equivDown(chPC, pc))
                continue;

            addSymmetrizationIndex(pc);

            // set PCs for channel's daughters
            for (auto& pcDaughPC : pc->daughters())
                for (const std::shared_ptr<Particle>& chDaugh : daughters())
                    if (std::dynamic_pointer_cast<DecayingParticle>(chDaugh))
                        for (auto& chDaughPC : std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->particleCombinations())
                            if (ParticleCombination::equivDown(pcDaughPC, chDaughPC))
                                std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->addSymmetrizationIndex(pcDaughPC);

        }
    }

    // next level
    for (auto d : daughters())
        d->setSymmetrizationIndexParents();
}

//--------------------------
DataAccessorSet DecayChannel::dataAccessors()
{
    DataAccessorSet V = {BlattWeisskopf_, SpinAmplitude_};

    for (auto& d : Daughters_) {
        // add daughter
        if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
            V.emplace(std::dynamic_pointer_cast<DecayingParticle>(d));
            // and its data accessors
            auto v = std::dynamic_pointer_cast<DecayingParticle>(d)->dataAccessors();
            V.insert(v.begin(), v.end());
        }
    }
    return V;
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
