#include "DecayingParticle.h"

#include "HelicitySpinAmplitude.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombinationCache.h"

#include <iomanip>
#include <memory>
#include <stdexcept>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    AmplitudeComponent(),
    Particle(q, mass, name),
    DataAccessor(),
    RadialSize_(std::make_shared<RealParameter>(radialSize)),
    Amplitude_(std::make_shared<ComplexCachedDataValue>(this))
{
}

//-------------------------
std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (Amplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        std::complex<double> a = Complex_0;

        /// \todo Is this the best way to do it? (loop over pc's then channels. or channels then pc's?)

        // sum up DecayChannel::amplitude over each channel
        for (auto& c : channels()) {
            if (c->hasSymmetrizationIndex(pc))
                a += c->amplitude(d, pc, dataPartitionIndex);
        }

        Amplitude_->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("DecayingParticle::amplitude - calculate amplitude for " << name() << " " << *pc << " = " << a);
        return a;
    }

    DEBUG("DecayingParticle::amplitude - using cached amplitude for " << name() << " " << *pc << " = " << Amplitude_->value(d, symIndex));
    return Amplitude_->value(d, symIndex);
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool C = DataAccessor::consistent();
    C &= Particle::consistent();

    if (RadialSize_->value() <= 0.) {
        FLOG(ERROR) << "Radial size not positive.";
        C &= false;
    }

    if (Channels_.empty()) {
        FLOG(ERROR) << "no channels specified.";
        return false; // further checks require at least one channel
    }

    // check no channel is empty
    if (std::any_of(Channels_.begin(), Channels_.end(), [](const std::unique_ptr<DecayChannel>& dc) {return !dc;})) {
        FLOG(ERROR) << "DecayChannel vector contains nullptr";
        C &= false;
    }
    // check all channels' parents point to this
    if (std::any_of(Channels_.begin(), Channels_.end(), [&](const std::unique_ptr<DecayChannel>& dc) {return dc and dc->decayingParticle() != this;})) {
        FLOG(ERROR) << "DecayChannel vector contains channel not pointing back to this";
        C &= false;
    }
    // check consistency of all channels
    std::for_each(Channels_.begin(), Channels_.end(), [&](const std::unique_ptr<DecayChannel>& dc) {if (dc) C &= dc->consistent();});

    // check if all channels lead to same final state particles
    /// \todo This isn't necessary, we should think how to change this. Example: D -> KKpipi, with f0->KK and f0->pipi
    std::vector<std::shared_ptr<FinalStateParticle> > fsps0 = finalStateParticles(0);
    std::sort(fsps0.begin(), fsps0.end());
    for (unsigned i = 1; i < nChannels(); ++i) {
        std::vector<std::shared_ptr<FinalStateParticle> > fsps = finalStateParticles(i);
        std::sort(fsps.begin(), fsps.end());
        if (fsps != fsps0) {
            LOG(ERROR) << "DecayingParticle::consistent() - final state of channel " << i << " does not match.";
            C &= false;
        }
    }

    return C;
}

//-------------------------
void DecayingParticle::addChannel(std::unique_ptr<DecayChannel> c)
{
    if (!c)
        throw exceptions::DecayChannelEmpty();

    if (c->particleCombinations().empty())
        throw exceptions::ParticleCombinationsEmpty();

    // check ISP
    if (initialStateParticle() and c->initialStateParticle() != initialStateParticle())
        throw exceptions::InitialStateParticleMismatch();

    Channels_.emplace_back(std::move(c));
    Channels_.back()->setDecayingParticle(this);

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations())
        addSymmetrizationIndex(pc);

    // add dependencies
    Amplitude_->addDependencies(Channels_.back()->ParametersItDependsOn());
    Amplitude_->addDependencies(Channels_.back()->CachedDataValuesItDependsOn());
}

//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (!Channels_.at(i))
        return std::vector<std::shared_ptr<FinalStateParticle>>();
    return Channels_[i]->finalStateParticles();
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
        }

        // set dependencies
        channel(i)->addSpinAmplitudeDependencies();

        // recurse down
        for (std::shared_ptr<Particle> d : channel(i)->daughters()) {
            if (std::dynamic_pointer_cast<DecayingParticle>(d))
                std::static_pointer_cast<DecayingParticle>(d)->optimizeSpinAmplitudeSharing();
        }
    }

    if (first) {
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
