#include "DecayChannel.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "Particle.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <assert.h>

namespace yap {

//-------------------------
DecayChannel::DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude) :
    DecayChannel( {daughterA, daughterB}, spinAmplitude)
{
}

//-------------------------
DecayChannel::DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude) :
    AmplitudeComponentDataAccessor(spinAmplitude->initialStateParticle()),
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
Amp DecayChannel::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    /// \todo check
    Amp a = BlattWeisskopf_.amplitude(d, pc) * SpinAmplitude_->amplitude(d, pc);
    for (auto& daugh : Daughters_) {
        a *= daugh->amplitude(d, pc);
    }

    DEBUG("DecayChannel: amplitude = " << a);

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
double DecayChannel::breakupMomentum() const
{
    if (Daughters_.size() != 2) {
        LOG(ERROR) << "DecayChannel::breakupMomentum() - channel has != 2 daughters. Cannot calculate!";
        return 0;
    }

    /// \todo take masses from mass shape instead?
    double m2_R =  pow(Parent_->mass(), 2);
    double m_a = Daughters_[0]->mass();
    double m_b = Daughters_[1]->mass();

    if (m_a == m_b) {
        return m2_R / 4.0 - m_a * m_a;
    }

    return (m2_R - (m_a + m_b) * (m_a + m_b)) *
           (m2_R - (m_a - m_b) * (m_a - m_b)) / m2_R / 4.0;
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
            std::vector<std::shared_ptr<FinalStateParticle> > ddaughters = std::dynamic_pointer_cast<DecayingParticle>(d)->finalStateParticles();
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

//-------------------------
void DecayChannel::clearSymmetrizationIndices()
{
    DataAccessor::clearSymmetrizationIndices();
    BlattWeisskopf_.clearSymmetrizationIndices();
    SpinAmplitude_->clearSymmetrizationIndices();
}

//-------------------------
void DecayChannel::setSymmetrizationIndexParents()
{
    //std::cout << "DecayChannel::setSymmetrizationIndexParents()\n";

    std::vector<std::shared_ptr<ParticleCombination> > chPCs = particleCombinations();

    // clean up PCs without parents
    std::vector<std::shared_ptr<ParticleCombination> > chPCsParents = particleCombinations();
    auto it = chPCsParents.begin();
    while (it != chPCsParents.end()) {
        if (not (*it)->parent()) {
            it = chPCsParents.erase(it);
        } else
            ++it;
    }
    clearSymmetrizationIndices();

    for (auto& pc : chPCsParents) {
        /*std::cout << "  add " << std::string(*pc);
        if  (pc->parent()){
          std::cout << " from decay " << std::string(*pc->parent());
        }
        std::cout  << " to channel " << std::string(*this) << "\n";*/

        addSymmetrizationIndex(pc);
    }


    for (auto& chPC : chPCs) {
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivDown(chPC, pc)) {

                //std::cout << std::string(*chPC) << " == " << std::string(*pc) << "\n";
                /*std::cout << "  add " << std::string(*pc);
                if  (pc->parent()){
                  std::cout << " from decay " << std::string(*pc->parent());
                }
                std::cout << " to channel " << std::string(*this) << "\n";*/

                addSymmetrizationIndex(pc);

                // set PCs for channel's daughters
                for (auto& pcDaughPC : pc->daughters()) {
                    for (const std::shared_ptr<Particle>& chDaugh : daughters()) {
                        if (std::dynamic_pointer_cast<DecayingParticle>(chDaugh))
                            for (auto& chDaughPC : std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->particleCombinations()) {
                                if (ParticleCombination::equivDown(pcDaughPC, chDaughPC)) {
                                    std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->addSymmetrizationIndex(pcDaughPC);
                                    /*std::cout << "   add " << std::string(*pcDaughPC);
                                    if (pcDaughPC->parent()) {
                                      std::cout << " from decay " << std::string(*pcDaughPC->parent());
                                    }
                                    std::cout << " to particle " << chDaugh->name() << "\n";*/
                                }
                            }
                    }
                }
            }
        }
    }

    /*std::cout << "Particle combinations in channel " <<  std::string(*this) << "\n";
    for (auto& chPC : particleCombinations()) {
      std::cout << "  " << std::string(*chPC);
      if (chPC->parent())
        std::cout << " from decay " << std::string(*chPC->parent());
      std::cout << "\n";
    }

    for (auto d : daughters()) {
      if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
        std::cout << "    PCs in particle " << d->name() << ":\n";
        for (auto& meh : std::dynamic_pointer_cast<DecayingParticle>(d)->particleCombinations()) {
            std::cout << "      " << std::string(*meh);
            if (meh->parent())
              std::cout << " from decay " << std::string(*meh->parent());
            std::cout <<"\n";
        }
      }
    }*/

    // next level
    for (auto d : daughters())
        d->setSymmetrizationIndexParents();

}


}
