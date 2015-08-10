#include "DecayingParticle.h"

#include "FinalStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"

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
        return false; // furhter checks require at least one channel
    }

    const std::vector<const FinalStateParticle*> fsps0 = this->finalStateParticles(0);

    for (unsigned i = 0; i < nChannels(); ++i) {
        const DecayChannel* c = this->getChannel(i);
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
    const Daughters& daughters = this->getChannel(channel)->daughters();

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

}
