#include "Resonance.h"

namespace yap {


//-------------------------
const std::vector<const FinalStateParticle*> Resonance::getFinalStateParticles(unsigned int channel) const {
  std::vector<const FinalStateParticle*> fsps;
  const Daughters& daughters = Channels_.at(channel)->getDaughters();

  for (const Particle* d : daughters) {
    if (dynamic_cast<const FinalStateParticle*>(d))
      fsps.push_back(static_cast<const FinalStateParticle*>(d));
    else {
      const std::vector<const FinalStateParticle*> ddaughters =
          dynamic_cast<const Resonance*>(d)->getFinalStateParticles();
      fsps.insert(fsps.end(), ddaughters.begin(), ddaughters.end());
    }
  }

  return fsps;
}

}
