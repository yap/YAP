#include "DecayChannel.h"

#include "container_utils.h"
#include "DecayingParticle.h"
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
DecayChannel::DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude, DecayingParticle* parent) :
    DataAccessor(),
    Parent_(parent),
    Daughters_(daughters),
    BlattWeisskopf_(std::make_unique<BlattWeiskopf>(this)),
    SpinAmplitude_(spinAmplitude),
    FreeAmplitude_(new ComplexParameter(1.)),
    FixedAmplitude_(new ComplexCachedDataValue(this))
{
    // check daughter size
    if (daughters.empty())
        throw std::runtime_error("daughters is empty");
    if (daughters.size() == 1)
        throw std::runtime_error("only one daughter provided");
    if (daughters.size() > 2)
        throw std::runtime_error("currently only supporting two-body decays");

    // check that first daughter's ISP is not nullptr
    if (daughters[0]->initialStateParticle() == nullptr)
        throw std::runtime_error("daughters' initial-state particle is unset");
    // and that all have same ISP (trivially checks 0th against itself)
    for (auto& d : daughters)
        if (d->initialStateParticle() != daughters[0]->initialStateParticle())
            throw std::runtime_error("daughters' initial-state particles not all the same.");

    // set intial-state particle (without adding this to the ISP yet)
    setInitialStateParticle(daughters[0]->initialStateParticle());
    BlattWeiskopf_->setInitialStateParticle(initialStateParticle());
    SpinAmplitude_->setInitialStateParticle(initialStateParticle());

    /// set dependencies
    FixedAmplitude_->addDependencies(BlattWeisskopf_->ParametersItDependsOn());
    FixedAmplitude_->addDependencies(BlattWeisskopf_->CachedDataValuesItDependsOn());

    // Spin amplitude dependencies are added via
    // addSpinAmplitudeDependencies() after sharing SpinAmplitudes

    // add daughter dependencies to FixedAmplitude_
    for (int i = 0; i < int(Daughters_.size()); ++i)
        if (auto d = std::dynamic_pointer_cast<DecayingParticle>(Daughters_[i]))
            for (auto& c : d->CachedDataValuesItDependsOn())
                FixedAmplitude_->addDependency(c, i);

    // set symmetrization indices
    std::vector<ParticleCombinationVector> PCs;
    for (std::shared_ptr<Particle> d : Daughters_)
        PCs.push_back(d->particleCombinations());

    /// \todo remove hardcoding for two daughters so applies to n daughters
    for (auto& PCA : PCs[0]) {
        for (auto& PCB : PCs[1]) {

            // check that PCA and PCB don't overlap in FSP content
            if (overlap(PCA->indices(), PCB->indices()))
                continue;

            // for identical particles, check if swapped particle combination is already added
            if (Daughters_[0] == Daughters_[1]) {
                // get (B,A) combination from cache
                auto b_a = intialStateParticle::particleCombinationCache.find(new ParticleCombination({PCB, PCA}));
                // if b_a is not in cache, it can't be in SymmetrizationIndices_
                if (b_a and SymmetrizationIndices_.find(b_a))
                    // if (B,A) already added, don't proceed to adding (A,B)
                    continue;
            }

            // add (A,B)
            addSymmetrizationIndex(initialStateParticle::particleCombinationCache[ {PCA, PCB}]);
        }
    }
}

//-------------------------
std::complex<double> DecayChannel::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    DEBUG("DecayChannel::amplitude - " << std::string(*this) << " " << *pc);

    unsigned symIndex = symmetrizationIndex(pc);

    if (FixedAmplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        std::complex<double> a = BlattWeisskopf_->amplitude(d, pc, dataPartitionIndex) * SpinAmplitude_->amplitude(d, pc, dataPartitionIndex);

        auto& pcDaughters = pc->daughters();
        for (unsigned i = 0; i < Daughters_.size(); ++i)
            a *= Daughters_[i]->amplitude(d, pcDaughters.at(i), dataPartitionIndex);

        FixedAmplitude_->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("DecayChannel::amplitude - calculated fixed amplitude for " << std::string(*this) << " " << *pc << " = " << a);
        return FreeAmplitude_->value() * a;
    }

    DEBUG("DecayChannel::amplitude - use cached fixed amplitude for " << std::string(*this) << " " << *pc << " = " << FixedAmplitude_->value(d, symIndex));
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
                    [&](std::shared_ptr<ParticleCombination> pc){return pc->daughters().size() != Daughters_.size();})) {
        FLOG(ERROR) << "DecayChannel and its particleCombinations do not have the same number of daughters.";
        C &= false;
    }

    // check no daughters is empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [](std::shared_ptr<Particle> d){return !d;})) {
        FLOG(ERROR) << "null pointer in daughters vector.";
        C &= false;
    }
    // check daughters
    std::for_each(Daughters_.begin(), Daughters_.end(), [&](std::shared_ptr<Particle> d){C &= (d) ? d->consistent() : false;});
    
    // Check Blatt-Weisskopf object
    if (!BlattWeiskopf) {
        FLOG(ERROR) << "BlattWeiskopf is nullptr";
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
        if (SpinAmplitude_->initialQuantumNumbers() != parent()->quantumNumbers()) {
            FLOG(ERROR) << "quantum numbers of parent " << parent()->quantumNumbers()
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
                                      [](double m, std::shared_ptr<Particle> d){return m + (d ? d->mass()->value() : 0);});
    if (mass_sum > parent()->mass()->value()) {
        FLOG(ERROR) << "sum of daughter's masses (" << mass_sum << ")"
                    << "is bigger than resonance mass (" << parent()->mass()->value() << ").";
        C &= false;
    }

    return C;
}

//-------------------------
std::string to_string(const DecayChannel& dc)
{
    std::string s;
    if (dc.parent())
        s += dc.parent()->name() + " ->";
    if (dc.daughters().empty())
        s += " (nothing)";
    else
        for (auto& d : dc.daughters())
            s += " " + d->name();
    if (dc.spinAmplitude())
        s += " " + std::string(*SpinAmplitude_);
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
            throw std::runtime_error("invalid daughter");
        }
    }

    return fsps;
}

//-------------------------
void DecayChannel::setInitialStateParticle(InitialStateParticle* isp)
{
    DataAccessor::setInitialStateParticle(isp);
    if (!SpinAmplitude_)
        std::throw std::runtime_error("SpinAmplitude not set.");
    SpinAmplitude_->setInitialStateParticle(isp);
    if (!BlattWeisskopf)
        std::throw std::runtime_error("BlattWeiskopf not set.");
    BlattWeisskopf_->setInitialStateParticle(isp);

    // hand ISP to daughters
    for (auto d : Daughters_)
        if (std::dynamic_pointer_cast<DecayingParticle>(d))
            std::static_pointer_cast<DecayingParticle>(d)->setInitialStateParticle(isp);
}

//-------------------------
void DecayChannel::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    DataAccessor::addSymmetrizationIndex(c);
    BlattWeisskopf_->addSymmetrizationIndex(c);
    SpinAmplitude_->addSymmetrizationIndex(c);
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

    for (auto& pc : chPCsParents) {
        addSymmetrizationIndex(pc);
    }


    for (auto& chPC : chPCs) {
        for (auto& wpc : ParticleCombination::cache) {
            if (wpc.expired())
                continue;
            auto pc = wpc.lock();
            if (ParticleCombination::equivDown(chPC, pc)) {

                addSymmetrizationIndex(pc);

                // set PCs for channel's daughters
                for (auto& pcDaughPC : pc->daughters()) {
                    for (const std::shared_ptr<Particle>& chDaugh : daughters()) {
                        if (std::dynamic_pointer_cast<DecayingParticle>(chDaugh))
                            for (auto& chDaughPC : std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->particleCombinations()) {
                                if (ParticleCombination::equivDown(pcDaughPC, chDaughPC)) {
                                    std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->addSymmetrizationIndex(pcDaughPC);
                                }
                            }
                    }
                }
            }
        }
    }

    // next level
    for (auto d : daughters())
        d->setSymmetrizationIndexParents();

}

//-------------------------
void DecayChannel::addSpinAmplitudeDependencies()
{
    FixedAmplitude_->addDependencies(SpinAmplitude_->ParametersItDependsOn());
    FixedAmplitude_->addDependencies(SpinAmplitude_->CachedDataValuesItDependsOn());
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
