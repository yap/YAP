#include "DecayChannel.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "Particle.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <assert.h>

namespace yap {

//-------------------------
DecayChannel::DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude, DecayingParticle* parent) :
    DecayChannel( {daughterA, daughterB}, spinAmplitude, parent)
{
}

//-------------------------
DecayChannel::DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude, DecayingParticle* parent) :
    AmplitudeComponentDataAccessor(),
    Parent_(parent),
    Daughters_(daughters),
    BlattWeisskopf_(nullptr), // see comment below!
    SpinAmplitude_(spinAmplitude),
    FreeAmplitude_(new ComplexParameter(0)),
    FixedAmplitude_(new CachedDataValue(this, {})),
    NominalBreakupMomentum2_(new RealCachedValue())
{
    /// this is done here because BlattWeisskopf needs a constructed DecayChannel object to set its dependencies
    std::unique_ptr<BlattWeisskopf> bw(new BlattWeisskopf(this));
    BlattWeisskopf_.swap(bw);

    // NominalBreakupMomentum2_ dependencies
    NominalBreakupMomentum2_->addDependency(Parent_->mass());
    for (auto& d : Daughters_)
        NominalBreakupMomentum2_->addDependency(d->mass());

    /// \todo set all other dependencies


    // set symmetrization indices
    std::vector<std::vector<std::shared_ptr<const ParticleCombination> > > PCs;
    for (std::shared_ptr<Particle> d : Daughters_) {
        if (std::dynamic_pointer_cast<DataAccessor>(d))
            PCs.push_back(std::dynamic_pointer_cast<DataAccessor>(d)->particleCombinations());
        else if (std::dynamic_pointer_cast<FinalStateParticle>(d))
            PCs.push_back(std::static_pointer_cast<FinalStateParticle>(d)->particleCombinations());
        else
            LOG(ERROR) << "DecayChannel() - cannot get ParticleCombinations from daughter " << d->name();
    }

    if (PCs.size() < 2) {
        LOG(ERROR) << "DecayChannel::DecayChannel - too few daughters provided.";
        return;
    }

    if (PCs.size() != 2) {
        LOG(ERROR) << "DecayChannel::DecayChannel - currently only accepting two-body decays.";
        return;
    }

    /// \todo how to for three?
    // hard-coded for two
    for (std::shared_ptr<const ParticleCombination> PCA : PCs[0])
        for (std::shared_ptr<const ParticleCombination> PCB : PCs[1])
            if (!PCA->sharesIndices(PCB)) {
                std::shared_ptr<const ParticleCombination> a_b = ParticleCombination::uniqueSharedPtr({PCA, PCB});

                bool can_has_symmetrization = true;
                if (Daughters_[0] == Daughters_[1]) {
                    ParticleCombination b_a = ParticleCombination({PCB, PCA});
                    for (auto& kv : SymmetrizationIndices_)
                        if (*(kv.first) == b_a) {
                            can_has_symmetrization = false;
                            break;
                        }
                }
                if (can_has_symmetrization)
                    addSymmetrizationIndex(a_b);
            }
}

//-------------------------
std::complex<double> DecayChannel::calcAmplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo check
    std::complex<double> a = BlattWeisskopf_->amplitude(d, pc) * SpinAmplitude_->amplitude(d, pc);

    if (a == Complex_0)
        return a;

    auto& pcDaughters = pc->daughters();
    for (unsigned i = 0; i < Daughters_.size(); ++i) {
        a *= Daughters_[i]->amplitude(d, pcDaughters.at(i));
    }

    DEBUG("DecayChannel " << std::string(*this) << " amplitude for " << std::string(*pc) << " = " << a);

    return a;
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool result = true;

    result &= AmplitudeComponentDataAccessor::consistent();
    if (!result) {
        LOG(ERROR) << "Channel's AmplitudeComponentDataAccessor is not consistent:  " << static_cast<std::string>(*this) << "\n";
    }

    // check number of daughters greater than 1
    if (Daughters_.size() < 2) {
        LOG(ERROR) << "DecayChannel::consistent() - invalid number of daughters (" << Daughters_.size() << " < 2).";
        result = false;
    }

    // currently only allowing exactly two daughters
    /// \todo allow more than two daugters?
    if (Daughters_.size() != 2) {
        LOG(ERROR) << "DecayChannel::consistent() - invalid number of daughters (" << Daughters_.size() << " != 2).";
        result = false;
    }

    // compare number of daughters
    for (auto& pc : particleCombinations())
        if (Daughters_.size() != pc->daughters().size()) {
            LOG(ERROR) << "DecayChannel::consistent() - DecayChannel and its particleCombinations do not have the same number of daughters.";
            result = false;
        }

    // check daughters
    bool prevResult = result;
    for (std::shared_ptr<Particle> d : Daughters_)  {
        if (!d) {
            LOG(ERROR) << "DecayChannel::consistent() - null pointer in daughters vector.";
            result = false;
        } else
            result &= d->consistent();
    }
    if (prevResult != result)
        LOG(ERROR) << "DecayChannel::consistent() - daughter(s) inconsistent";

    // Check Blatt-Weisskopf object
    result &= BlattWeisskopf_->consistent();
    // check if BlattWeisskopf points back to this DecayChannel
    if (this != BlattWeisskopf_->decayChannel()) {
        LOG(ERROR) << "DecayChannel::consistent() - BlattWeisskopf does not point back to this DecayChannel.";
        result =  false;
    }


    // Check SpinAmplitude object
    if (!SpinAmplitude_) {
        LOG(ERROR) << "DecayChannel::consistent() - no SpinAmplitude object set.";
        result = false;
    } else {
        result &= SpinAmplitude_->consistent();

        // check size of spin amplitude quantum numbers and size of daughters
        if (SpinAmplitude_->finalQuantumNumbers().size() != Daughters_.size()) {
            LOG(ERROR) << "DecayChannel::consistent() - quantum numbers object and daughters object size mismatch";
            result = false;
        }

        // check if QuantumNumbers of SpinAmplitude objects match with Particles
        if (SpinAmplitude_->initialQuantumNumbers() != parent()->quantumNumbers()) {
            LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of parent "
                       << parent()->quantumNumbers() << " and SpinAmplitude "
                       << SpinAmplitude_->initialQuantumNumbers() << " don't match.";
            result = false;
        }

        for (unsigned i = 0; i < Daughters_.size(); ++i) {
            if (SpinAmplitude_->finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
                LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of daughter " << i << " "
                           << Daughters_[i]->quantumNumbers() << " and SpinAmplitude "
                           << SpinAmplitude_->finalQuantumNumbers()[i] << " don't match.";
                result = false;
            }
        }
    }

    // check masses
    double finalMass = 0;
    for (std::shared_ptr<Particle> d : Daughters_)
        finalMass += (!d) ? 0 : d->mass()->value();
    if (finalMass > parent()->mass()->value()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        result =  false;
    }

    if (!result)
        LOG(ERROR) << "Channel is not consistent:  " << static_cast<std::string>(*this) << "\n";

    return result;
}

//-------------------------
std::shared_ptr<RealCachedValue> DecayChannel::breakupMomentum2()
{
    if (NominalBreakupMomentum2_->calculationStatus() == kUncalculated) {
        if (Daughters_.size() != 2) {
            LOG(FATAL) << "DecayChannel::breakupMomentum() - wrong number of daughters (" << Daughters_.size() << " != 2). Cannot calculate!";
            NominalBreakupMomentum2_->setValue(std::numeric_limits<double>::quiet_NaN());
        }

        /// \todo take masses from mass shape instead?
        /// or do we need the invariant masses instead of the nominal masses?
        double m2_R =  pow(Parent_->mass()->value(), 2);
        double m_a = Daughters_[0]->mass()->value();
        double m_b = Daughters_[1]->mass()->value();

        if (m_a == m_b) {
            NominalBreakupMomentum2_->setValue( m2_R / 4. - m_a * m_a );
        }

        NominalBreakupMomentum2_->setValue( (m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4. );
    }

    return NominalBreakupMomentum2_;
}

//-------------------------
DecayChannel::operator std::string() const
{
    std::string result;
    if (Parent_)
        result += Parent_->name() + " ->";
    if (Daughters_.empty())
        result += " (nothing)";
    for (std::shared_ptr<Particle> d : Daughters_)
        result += " " + d->name();
    if (SpinAmplitude_)
        result += " " + std::string(*SpinAmplitude_);
    return result;
}

//-------------------------
std::vector<std::shared_ptr<FinalStateParticle> > DecayChannel::finalStateParticles() const
{
    std::vector<std::shared_ptr<FinalStateParticle> > fsps;

    for (std::shared_ptr<Particle> d : Daughters_) {

        if (std::dynamic_pointer_cast<FinalStateParticle>(d)) {
            fsps.push_back(std::static_pointer_cast<FinalStateParticle>(d));

        } else if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
            std::vector<std::shared_ptr<FinalStateParticle> > ddaughters = std::dynamic_pointer_cast<DecayingParticle>(d)->finalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());

        } else {
            LOG(ERROR) << "DecayingParticle::finalStateParticles() - Daughter is neither a FinalStateParticle nor a DecayingParticle. DecayChannel is inconsistent.";
        }
    }

    return fsps;
}

//-------------------------
void DecayChannel::setInitialStateParticle(InitialStateParticle* isp)
{
    AmplitudeComponentDataAccessor::setInitialStateParticle(isp);
    SpinAmplitude_->setInitialStateParticle(initialStateParticle());
    BlattWeisskopf_->setInitialStateParticle(initialStateParticle());

    // hand ISP to daughters
    for (auto d : Daughters_)
        if (std::dynamic_pointer_cast<DecayingParticle>(d))
            std::static_pointer_cast<DecayingParticle>(d)->setInitialStateParticle(initialStateParticle());
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
    auto chPCs = particleCombinations();

    // clean up PCs without parents
    auto chPCsParents = particleCombinations();
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
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
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
void DecayChannel::precalculate()
{
    // calculates it
    breakupMomentum2();

    /// \todo find a solution to do this not recursively?
    BlattWeisskopf_->precalculate();
    SpinAmplitude_->precalculate();
    for (auto& d : Daughters_)
        d->precalculate();
}

//-------------------------
void DecayChannel::finishedPrecalculation()
{
    BlattWeisskopf_->finishedPrecalculation();
    SpinAmplitude_->finishedPrecalculation();
    for (auto& d : Daughters_)
        d->finishedPrecalculation();
}


}
