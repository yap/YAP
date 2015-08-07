#include "ParticleFactory.h"
#include "Resonance.h"

#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP

int main( int argc, char** argv)
{
  yap::ParticleFactory factory;

  // final state particles
  yap::Particle* piPlus1 = factory.createFinalStateParticle(211);
  yap::Particle* piPlus2 = factory.createFinalStateParticle(211);
  yap::Particle* piMinus1 = factory.createFinalStateParticle(-211);
  yap::Particle* piMinus2 = factory.createFinalStateParticle(-211);


}
