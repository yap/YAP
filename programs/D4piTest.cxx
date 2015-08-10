#include "ParticleFactory.h"
//#include "SpinAmplitude.h"
//#include "Resonance.h"

#include <assert.h>
#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP

int main( int argc, char** argv)
{
  yap::ParticleFactory factory;

  // final state particles
  yap::FinalStateParticle* piPlus = factory.createFinalStateParticle(211);
  yap::FinalStateParticle* piMinus = factory.createFinalStateParticle(-211);

  double radialSize = 1.;
  unsigned L = 1;
  yap::Resonance* rho = factory.createResonanceBreitWigner(113, radialSize);
  yap::SpinAmplitude* rhoAmplitude = new yap::SpinAmplitude(rho->quantumNumbers(), piPlus->quantumNumbers(), piMinus->quantumNumbers());
  rho->addChannel(yap::DecayChannel(piPlus, piMinus, L, *rhoAmplitude));

  assert(rho->consistent());
}
