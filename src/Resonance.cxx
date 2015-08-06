#include "Resonance.h"
#include "logging.h"

namespace yap {

//-------------------------
bool Resonance::consistent() const
{
    bool consistent = true;
    consistent &= MassShape_.consistent();

    if (RadialSize_ <= 0.) {
        LOG(ERROR) << "Resonance::consistent() - Radial size not positive.";
        consistent = false;
    }

    if (Channels_.empty()) {
        LOG(ERROR) << "Resonance::consistent() - no channels specified.";
        return false; // furhter checks require at least one channel
    }

    const std::vector<const FinalStateParticle*> fsps0 = this->finalStateParticles(0);

    for (DecayChannel c : Channels_) {
        if (this != c.resonance()) {
            LOG(ERROR) << "Resonance::consistent() - DecayChannels does not point back to this Resonance.";
            return false; // channel consistency check requires correct pointer
        }

        consistent &= c.consistent();
    }

    // check if all channels lead to same final state particles
    if (nChannels() > 1) {
        for (unsigned int i = 1; i < nChannels(); ++i) {
            const std::vector<const FinalStateParticle*> fsps = this->finalStateParticles(0);
            if (fsps0.size() != fsps.size()) {
                LOG(ERROR) << "Resonance::consistent() - number of final state particles of different channels do not match.";
                consistent = false;
            }
            for (unsigned int j = 0; j < fsps0.size(); ++j) {
                // compare adresses to check if final state particles are really the same objects
                if (fsps0[j] != fsps[j]) {
                    LOG(ERROR) << "Resonance::consistent() - final state particles of different channels are not the same (objects).";
                    consistent = false;
                }
            }
        }
    }


    return consistent;
}

//-------------------------
const std::vector<const FinalStateParticle*> Resonance::finalStateParticles(unsigned int channel) const
{
    std::vector<const FinalStateParticle*> fsps;
    const Daughters& daughters = Channels_.at(channel).daughters();

    for (const Particle* d : daughters) {
        if (dynamic_cast<const FinalStateParticle*>(d))
            fsps.push_back(static_cast<const FinalStateParticle*>(d));
        else {
            const std::vector<const FinalStateParticle*> ddaughters =
                dynamic_cast<const Resonance*>(d)->finalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());
        }
    }

    return fsps;
}

}
