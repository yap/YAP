#include "DecayingParticle.h"

#include "CanonicalSpinAmplitude.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"

#include <iomanip>
#include <memory>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    Particle(q, mass, name),
    DataAccessor(),
    RadialSize_(radialSize)
{}

//-------------------------
Amp DecayingParticle::amplitude(DataPoint& d)
{
    // \todo implement
    return Amp(1);
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool result = true;

    result &= DataAccessor::consistent();
    result &= Particle::consistent();

    if (RadialSize_ <= 0.) {
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
void DecayingParticle::addChannel(DecayChannel* c)
{
    Channels_.push_back(std::unique_ptr<yap::DecayChannel>(c));
    Channels_.back()->setParent(this);

    if (c->particleCombinations().empty())
        LOG(ERROR) << "DecayingParticle::addChannel(c) - c->particleCombinations().empty()";

    for (std::shared_ptr<ParticleCombination> pc : c->particleCombinations())
        this->addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(pc));
}

//-------------------------
void DecayingParticle::addChannel(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned twoL)
{
    addChannel(new DecayChannel(A, B, std::make_shared<CanonicalSpinAmplitude>(quantumNumbers(), A->quantumNumbers(), B->quantumNumbers(), twoL)));
}

//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (!Channels_[i])
        return std::vector<std::shared_ptr<FinalStateParticle>>();
    return Channels_[i]->finalStateParticles();
}

//-------------------------
void DecayingParticle::optimizeSpinAmplitudeSharing()
{
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;
    static SharedSpinAmplitudeComparator comp;

    for (unsigned i = 0; i < nChannels(); ++i) {

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();

        bool found(false);
        for (auto sa : ampSet) {
            if (comp(sa, amp)) {
                LOG(INFO) << "Sharing spin amplitude " << std::string(*amp) << " and " << std::string(*sa);
                channel(i)->spinAmplitude() = sa;
                found = true;
                break;
            }
        }
        if (!found)
            ampSet.insert(amp);

        for (std::shared_ptr<Particle> d : channel(i)->daughters()) {
            if (std::dynamic_pointer_cast<DecayingParticle>(d))
                std::static_pointer_cast<DecayingParticle>(d)->optimizeSpinAmplitudeSharing();
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
            std::cout << "\n" << std::setw(level * (padding * 3 + 7 + paddingSpinAmp)) << "";

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
    static size_t paddingSpinAmp = 0;
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
            std::cout << "\n" << std::setw(level * (padding * 3 + 7 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();
        if (ampSet.find(amp) == ampSet.end()) {
            ampSet.insert(amp);
            std::cout << std::left << std::setw(paddingSpinAmp) << std::string(*amp);
        } else {
            std::cout << std::left << std::setw(paddingSpinAmp - 1) << "sharedSpinAmp";
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


}
