#include "ParticleFactory.h"
//#include "SpinAmplitude.h"
//#include "Resonance.h"

#include <assert.h>
#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP

int main( int argc, char** argv)
{
  yap::ParticleFactory factory("evt.pdl");

  // final state particles
  yap::FinalStateParticle* piPlus = factory.createFinalStateParticle(211);
  yap::FinalStateParticle* piMinus = factory.createFinalStateParticle(-211);

  // initial state particle
  double radialSize = 1.;
  yap::InitialStateParticle* D = factory.createInitialStateParticle(421, radialSize);

  // a_1 channels
  yap::Resonance* sigma = factory.createResonanceBreitWigner(9000221, radialSize);

  // rho rho
  yap::Resonance* rho = factory.createResonanceBreitWigner(113, radialSize);
  factory.createChannel(rho, piPlus, piMinus, 1);

  factory.createChannel(D, rho, rho, 0);
  factory.createChannel(D, rho, rho, 1);
  factory.createChannel(D, rho, rho, 2);

  // omega omega
  yap::Resonance* omega = factory.createResonanceBreitWigner(223, radialSize);
  factory.createChannel(omega, piPlus, piMinus, 1);

  factory.createChannel(D, omega, omega, 0);
  factory.createChannel(D, omega, omega, 1);
  factory.createChannel(D, omega, omega, 2);

  // rho omega
  factory.createChannel(D, rho, omega, 0);
  factory.createChannel(D, rho, omega, 1);
  factory.createChannel(D, rho, omega, 2);





  //std::cout << "rho " << rho->quantumNumbers();
  assert(D->consistent());
  std::cout << "alright! \n";
}
