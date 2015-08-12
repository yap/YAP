#include "DecayingParticle.h"

#include "FinalStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"

#include <iomanip>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    Particle(q, mass, name),
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
    bool consistent = true;

    consistent &= Particle::consistent();

    if (RadialSize_ <= 0.) {
        LOG(ERROR) << "DecayingParticle::consistent() - Radial size not positive.";
        consistent = false;
    }

    if (Channels_.empty()) {
        LOG(ERROR) << "DecayingParticle::consistent() - no channels specified.";
        return false; // further checks require at least one channel
    }

    const std::vector<const FinalStateParticle*> fsps0 = this->finalStateParticles(0);

    for (unsigned i = 0; i < nChannels(); ++i) {
        const DecayChannel* c = this->channel(i);
        if (this != c->parent()) {
            LOG(ERROR) << "DecayingParticle::consistent() - DecayChannels does not point back to this DecayingParticle.";
            return false; // channel consistency check requires correct pointer
        }

        consistent &= c->consistent();
    }

    // check if all channels lead to same final state particles
    if (nChannels() > 1) {
        for (unsigned int i = 1; i < nChannels(); ++i) {
            const std::vector<const FinalStateParticle*> fsps = this->finalStateParticles(0);
            if (fsps0.size() != fsps.size()) {
                LOG(ERROR) << "DecayingParticle::consistent() - number of final state particles of different channels do not match.";
                consistent = false;
            }
            for (unsigned int j = 0; j < fsps0.size(); ++j) {
                // compare adresses to check if final state particles are really the same objects
                if (fsps0[j] != fsps[j]) {
                    LOG(ERROR) << "DecayingParticle::consistent() - final state particles of different channels are not the same (objects).";
                    consistent = false;
                }
            }
        }
    }


    return consistent;
}

//-------------------------
const std::vector<const FinalStateParticle*> DecayingParticle::finalStateParticles(unsigned int channel) const
{
    std::vector<const FinalStateParticle*> fsps;
    const Daughters& daughters = this->channel(channel)->daughters();

    for (const Particle* d : daughters) {
        if (dynamic_cast<const FinalStateParticle*>(d))
            fsps.push_back(static_cast<const FinalStateParticle*>(d));
        else if (dynamic_cast<const DecayingParticle*>(d)) {
            const std::vector<const FinalStateParticle*> ddaughters =
                static_cast<const DecayingParticle*>(d)->finalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());
        } else {
            LOG(ERROR) << "DecayingParticle::finalStateParticles() - Daughter is neither a FinalStateParticle nor a DecayingParticle. DecayChannel is inconsistent.";
        }
    }

    return fsps;
}

//-------------------------
void DecayingParticle::addChannel(DecayChannel* c)
{
    Channels_.push_back(std::unique_ptr<yap::DecayChannel>(c));
    Channels_.back()->setParent(this);
}

//-------------------------
void DecayingParticle::optimizeSpinAmplitudeSharing()
{
    for (unsigned int i = 0; i < nChannels(); ++i) {
        // recursive call for daughters
        if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[0])) {
            static_cast<DecayingParticle*>(channel(i)->daughters()[0])->optimizeSpinAmplitudeSharing();
        }
        if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[1])) {
            static_cast<DecayingParticle*>(channel(i)->daughters()[1])->optimizeSpinAmplitudeSharing();
        }

        for (unsigned int j = i + 1; j < nChannels(); ++j) {
            // compare pointer adresses
            if ( channel(i)->spinAmplitude() == channel(j)->spinAmplitude() )
                continue; // same objects already

            // compare SpinAmplitude objects
            if  ( *(channel(i)->spinAmplitude()) == *(channel(j)->spinAmplitude()) ) {
                LOG(INFO) << "Share amplitudes of  "
                          << this->name() << " -> "
                          << channel(i)->daughters()[0]->name()
                          << " " << channel(i)->daughters()[1]->name()
                          << " (l=" << (int)channel(i)->spinAmplitude()->decayAngularMomentum() << ")"
                          << "  and  " << this->name() << " -> "
                          << channel(j)->daughters()[0]->name()
                          << " " << channel(j)->daughters()[1]->name()
                          << " (l=" << (int)channel(j)->spinAmplitude()->decayAngularMomentum() << ")";

                channel(j)->sharedSpinAmplitude().reset();
                channel(j)->sharedSpinAmplitude() = channel(i)->sharedSpinAmplitude() ;
            }

        }
    }
}

//-------------------------
void DecayingParticle::printDecayChainLevel(int level) const
{
    // get maximum length of particle names
    static unsigned padding = 0;
    if (padding == 0 || level == -1) {
        for (unsigned int i = 0; i < nChannels(); ++i) {
            padding = std::max(padding, (unsigned)this->name().length());
            padding = std::max(padding, (unsigned)channel(i)->daughters()[0]->name().length());
            padding = std::max(padding, (unsigned)channel(i)->daughters()[1]->name().length());

            if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[0])) {
                static_cast<DecayingParticle*>(channel(i)->daughters()[0])->printDecayChainLevel(-1);
            }
            if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[1])) {
                static_cast<DecayingParticle*>(channel(i)->daughters()[1])->printDecayChainLevel(-1);
            }

            if (level == -1)
                return;
        }
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 13)) << "";
        std::cout << std::left << std::setw(padding) << this->name() << " -> "
                  << std::setw(padding) << channel(i)->daughters()[0]->name()
                  << " " << std::setw(padding) << channel(i)->daughters()[1]->name()
                  << "(l=" << (int)channel(i)->spinAmplitude()->decayAngularMomentum() << ")";
        if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[0])) {
            std::cout << ",  ";
            static_cast<DecayingParticle*>(channel(i)->daughters()[0])->printDecayChainLevel(level + 1);
        }
        if (dynamic_cast<DecayingParticle*>(channel(i)->daughters()[1])) {
            std::cout << ",  ";
            static_cast<DecayingParticle*>(channel(i)->daughters()[1])->printDecayChainLevel(level + 1);
        }
    }

    if (level == 0)
        std::cout << "\n";
}


}
