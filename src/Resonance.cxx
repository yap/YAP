#include "Resonance.h"

namespace yap {


//-------------------------
const std::vector<FinalStateParticle&> Resonance::getFinalStateParticles(unsigned int channel) const
{
    std::vector<FinalStateParticle&> fsps;
    const std::array<Particle&, 2>& daughters = channels_.at(channel).getDaughters();

    for (const Particle& d : daughters) {
        if (dynamic_cast<FinalStateParticle&>(d))
            fsps.push_back(d);
        else {
            const std::vector<FinalStateParticle&> ddaughters =
                dynamic_cast<Resonance&>(d).getFinalStateParticles();
            fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());
        }
    }

    return fsps;
}

}
