#include "DecayingParticle.h"

#include "HelicitySpinAmplitude.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"
#include "SpinUtilities.h"

#include <iomanip>
#include <memory>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    AmplitudeComponent(),
    Particle(q, mass, name),
    DataAccessor(),
    RadialSize_(new RealParameter(radialSize))
{}

//-------------------------
std::complex<double> DecayingParticle::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    // \todo implement and check
    /*std::complex<double> a = Complex_0;

    for (auto& c : channels()) {
        if (c->hasSymmetrizationIndex(pc))
            a += c->freeAmplitude()->value() * c->amplitude(d, pc);
    }

    DEBUG("DecayingParticle " << name() << ": amplitude for " << std::string(*pc) << " = " << a);

    return a;*/
    return Complex_0;
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool result = true;

    result &= DataAccessor::consistent();
    result &= Particle::consistent();

    if (RadialSize_->value() <= 0.) {
        LOG(ERROR) << "DecayingParticle::consistent() - Radial size not positive.";
        result = false;
    }

    if (Channels_.empty()) {
        LOG(ERROR) << "DecayingParticle::consistent() - no channels specified.";
        return false; // further checks require at least one channel
    }

    for (auto& c : Channels_) {
        if (!c) {
            LOG(ERROR) << "DecayingParticle::consistent() - DecayChannels contains null pointer.";
            result = false;
        } else if (this != c->parent()) {
            LOG(ERROR) << "DecayingParticle::consistent() - DecayChannels does not point back to this DecayingParticle.";
            result = false; // channel consistency check requires correct pointer
        }
        result &= c->consistent();
    }

    // check if all channels lead to same final state particles
    std::vector<std::shared_ptr<FinalStateParticle> > fsps0 = finalStateParticles(0);
    for (unsigned i = 1; i < nChannels(); ++i)
        if (finalStateParticles(i) != fsps0) {
            LOG(ERROR) << "DecayingParticle::consistent() - final state of channel " << i << " does not match.";
            result = false;
        }

    return result;
}

//-------------------------
void DecayingParticle::addChannel(std::unique_ptr<DecayChannel>& c)
{
    Channels_.push_back(std::unique_ptr<yap::DecayChannel>(nullptr));
    Channels_.back().swap(c);
    Channels_.back()->setInitialStateParticle(initialStateParticle());

    if (Channels_.back()->particleCombinations().empty())
        LOG(ERROR) << "DecayingParticle::addChannel(c) - c->particleCombinations().empty()";

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations())
        addSymmetrizationIndex(pc);
}

//-------------------------
void DecayingParticle::addChannels(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned maxTwoL)
{
    // loop over possible l
    for (unsigned twoL = 0; twoL < maxTwoL; ++twoL) {
        if (!SpinAmplitude::angularMomentumConserved(quantumNumbers(), A->quantumNumbers(), B->quantumNumbers(), twoL))
            continue;

        std::unique_ptr<DecayChannel> chan( new DecayChannel(A, B, std::make_shared<HelicitySpinAmplitude>(quantumNumbers(), A->quantumNumbers(), B->quantumNumbers(), twoL), this) );

        bool notZero(false);
        std::vector<std::shared_ptr<const ParticleCombination>> PCs;

        for (char twoLambda = -quantumNumbers().twoJ(); twoLambda <= quantumNumbers().twoJ(); twoLambda += 2) {
            for (auto& pc : chan->particleCombinations()) {
                std::shared_ptr<ParticleCombination> pcHel(new ParticleCombination(*pc));
                pcHel -> setTwoLambda(twoLambda);

                if (std::static_pointer_cast<HelicitySpinAmplitude>(chan->spinAmplitude())->calculateClebschGordanCoefficient(pcHel) != 0) {
                    PCs.push_back(ParticleCombination::uniqueSharedPtr(pcHel));
                    notZero = true;

                }
            }
        }

        if (notZero) {
            DEBUG("add channel " << std::string(*chan) << " to " << name());
            chan->clearSymmetrizationIndices();
            for (auto& pc : PCs) {
                chan->addSymmetrizationIndex(pc);
            }

            addChannel(chan);
        }
    }

}


//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (!Channels_.at(i))
        return std::vector<std::shared_ptr<FinalStateParticle>>();
    return Channels_[i]->finalStateParticles();
}

//-------------------------
void DecayingParticle::setInitialStateParticle(InitialStateParticle* isp)
{
    DataAccessor::setInitialStateParticle(isp);
    // hand ISP to channels
    for (auto& c : Channels_)
        c->setInitialStateParticle(initialStateParticle());
}

//-------------------------
void DecayingParticle::optimizeSpinAmplitudeSharing()
{
    static bool firstCalled(true);
    bool first = firstCalled;

    firstCalled = false;

    /// \todo Will not work if we have more than one initial state particle
    /// make ampSet member of initialStateParticle ???
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;
    static SharedSpinAmplitudeComparator comp;

    for (unsigned i = 0; i < nChannels(); ++i) {

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();

        bool found(false);
        for (auto sa : ampSet) {
            if (comp(sa, amp)) {
                LOG(INFO) << "Sharing spin amplitude " << std::string(*sa) << " and " << std::string(*amp) << " (" << amp.get() << ")";
                channel(i)->spinAmplitude() = sa;
                found = true;
                break;
            }
        }
        if (!found) {
            ampSet.insert(amp);
            //LOG(INFO) << "Size of ampSet: " << ampSet.size();
        }

        for (std::shared_ptr<Particle> d : channel(i)->daughters()) {
            if (std::dynamic_pointer_cast<DecayingParticle>(d))
                std::static_pointer_cast<DecayingParticle>(d)->optimizeSpinAmplitudeSharing();
        }
    }

    if (first) {
        /*DEBUG("AmpSet:");
        for (auto& amp : ampSet) {
            DEBUG("  " << amp.get());
        }*/

        std::set<DataAccessor*> removeAmps;
        for (DataAccessor* dataAcc : initialStateParticle()->DataAccessors_) {
            if (dynamic_cast<SpinAmplitude*>(dataAcc)) {
                bool found(false);
                for (auto& amp : ampSet) {
                    if (amp.get() == dataAcc) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    removeAmps.insert(dataAcc);
                }
            }
        }

        for (DataAccessor* dataAcc : removeAmps) {
            initialStateParticle()->removeDataAccessor(dataAcc);
            LOG(INFO) << "remove unused SpinAmplitude " << dataAcc << " from InitialStateParticle's DataAccessors.";
        }
    }
}

//-------------------------
void DecayingParticle::printDecayChainLevel(int level) const
{
    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 0;
    if (padding == 0 || level == -1) {
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
            }
        }
        if (level == -1)
            return;
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 8 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();
        std::cout << std::left << std::setw(paddingSpinAmp)
                  << std::string(*(channel(i)->spinAmplitude()));

        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
}

//-------------------------
void DecayingParticle::printSpinAmplitudes(int level)
{
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;

    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 6;
    if (padding == 0 || level == -1) {
        ampSet.clear();
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
            }
        }
        if (level == -1)
            return;
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 8 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();
        if (ampSet.find(amp) == ampSet.end()) {
            ampSet.insert(amp);
            std::cout << std::left << std::setw(paddingSpinAmp) << std::string(*amp);
        } else {
            std::cout << std::left << std::setw(paddingSpinAmp - 1) << "shared";
        }


        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
}

//-------------------------
void DecayingParticle::setSymmetrizationIndexParents()
{
    //std::cout << "DecayingParticle::setSymmetrizationIndexParents()\n";

    // clean up PCs without parents
    auto PCsParents = particleCombinations();
    auto it = PCsParents.begin();
    while (it != PCsParents.end()) {
        if (!(*it)->parent()) {
            it = PCsParents.erase(it);
        } else
            ++it;
    }
    clearSymmetrizationIndices();

    for (auto& pc : PCsParents)
        addSymmetrizationIndex(pc);

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

}


}
