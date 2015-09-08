#include "DecayChannel.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "Particle.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude) :
    DecayChannel( {daughterA, daughterB}, spinAmplitude)
{
}

//-------------------------
DecayChannel::DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude) :
    AmplitudeComponent(),
    DataAccessor(),
    Daughters_(daughters),
    BlattWeisskopf_(this),
    SpinAmplitude_(spinAmplitude),
    FreeAmplitude_(0)
{
    // set symmetrization indices
    std::vector<std::vector<std::shared_ptr<ParticleCombination> > > PCs;
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
    for (std::shared_ptr<ParticleCombination> PCA : PCs[0])
        for (std::shared_ptr<ParticleCombination> PCB : PCs[1])
            if (!PCA->sharesIndices(PCB)) {
                std::shared_ptr<ParticleCombination> a_b = ParticleCombination::uniqueSharedPtr({PCA, PCB});
                //std::shared_ptr<ParticleCombination> a_b = std::make_shared<ParticleCombination>(ParticleCombination({PCA, PCB}));

                bool can_has_symmetrization = true;
                if (Daughters_[0] == Daughters_[1]) {
                    ParticleCombination b_a = ParticleCombination({PCB, PCA});
                    for (auto& kv : SymmetrizationIndices_)
                        if (*(kv.first) == b_a) {
                            can_has_symmetrization = false;
                            break;
                        }
                }
                if (can_has_symmetrization) {

                    /*bool newPC(false);
                    ParticleCombination* dummy = new ParticleCombination;
                    // set parents of daughter PCs
                    if (PCA->parent() == nullptr) {
                      // daughter has no parent yet. Set this as parent and add
                      PCA->setParent(a_b.get());
                    } else if (PCA->parent() == a_b.get()) {
                      // fine
                    } else {
                      // daughter has already different parent -> make copy, set this as parent and get unique shared_ptr
                      std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*PCA));
                      PCA.swap(copy);
                      PCA->setParent(dummy);
                      if (std::dynamic_pointer_cast<DataAccessor>(Daughters_[0]))
                          std::dynamic_pointer_cast<DataAccessor>(Daughters_[0])->addSymmetrizationIndex(PCA);
                      newPC = true;
                    }

                    if (PCB->parent() == nullptr) {
                      // daughter has no parent yet. Set this as parent and add
                      PCB->setParent(a_b.get());
                    } else if (PCB->parent() == a_b.get()) {
                      // fine
                    } else {
                      // daughter has already different parent -> make copy, set this as parent and get unique shared_ptr
                      std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*PCB));
                      PCB.swap(copy);
                      PCB->setParent(dummy);
                      if (std::dynamic_pointer_cast<DataAccessor>(Daughters_[1]))
                          std::dynamic_pointer_cast<DataAccessor>(Daughters_[1])->addSymmetrizationIndex(PCB);
                      newPC = true;
                    }

                    if (newPC) {
                      a_b = ParticleCombination::uniqueSharedPtr({PCA, PCB});
                      PCA->setParent(a_b.get());
                      PCB->setParent(a_b.get());
                    }*/

                    addSymmetrizationIndex(a_b);
                }
            }
}

//-------------------------
Amp DecayChannel::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    return Amp(1);
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool result = true;

    result &= DataAccessor::consistent();

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

    // check daughters
    for (std::shared_ptr<Particle> d : Daughters_)  {
        if (!d) {
            LOG(ERROR) << "DecayChannel::consistent() - null pointer in daughters vector.";
            result = false;
        } else
            result &= d->consistent();
    }
    if (!result)
        LOG(ERROR) << "DecayChannel::consistent() - daughter(s) inconsistent";

    // Check Blatt-Weisskopf object
    result &= BlattWeisskopf_.consistent();
    // check if BlattWeisskopf points back to this DecayChannel
    if (this != BlattWeisskopf_.decayChannel()) {
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
        finalMass += (!d) ? 0 : d->mass();
    if (finalMass > parent()->mass()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        result =  false;
    }

    if (!result)
        LOG(ERROR) << "Channel is not consistent:  " << static_cast<std::string>(*this) << "\n";

    return result;
}

//-------------------------
DecayChannel::operator std::string() const
{
    std::string result = parent()->name() + " ->";
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
            std::vector<std::shared_ptr<FinalStateParticle> > ddaughters = std::static_pointer_cast<DecayingParticle>(d)->finalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());

        } else {
            LOG(ERROR) << "DecayingParticle::finalStateParticles() - Daughter is neither a FinalStateParticle nor a DecayingParticle. DecayChannel is inconsistent.";
        }
    }

    return fsps;
}

//-------------------------
void DecayChannel::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
    DataAccessor::addSymmetrizationIndex(c);
    BlattWeisskopf_.addSymmetrizationIndex(c);
    SpinAmplitude_->addSymmetrizationIndex(c);
}


}
