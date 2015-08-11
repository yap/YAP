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
  yap::Resonance* rho = factory.createResonanceBreitWigner(113, radialSize);
  factory.createChannel(rho, piPlus, piMinus, 1);

  yap::InitialStateParticle* D = factory.createInitialStateParticle(421, radialSize);
  factory.createChannel(D, rho, rho, 0);
  factory.createChannel(D, rho, rho, 1);
  factory.createChannel(D, rho, rho, 2);



  assert(rho->consistent());
}
