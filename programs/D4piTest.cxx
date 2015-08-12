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


  // a_1 channels
  yap::Resonance* sigma = factory.createResonanceBreitWigner(9000221, radialSize);
  factory.createChannel(sigma, piPlus, piMinus, 0);

  yap::Resonance* a_1 = factory.createResonanceBreitWigner(20213, radialSize);
  factory.createChannel(a_1, sigma, piPlus, 1);

  factory.createChannel(a_1, rho, piPlus, 0); // S-wave
  factory.createChannel(a_1, rho, piPlus, 1); // not in Focus model
  factory.createChannel(a_1, rho, piPlus, 2); // D-wave

  factory.createChannel(D, a_1, piMinus, 1);


  // R pi pi channels
  //yap::Resonance* f_0_980 = factory.createResonanceBreitWigner(9000221, radialSize);
  //factory.createChannel(f_0_980, piPlus, piMinus, 0);



  /*std::cout << piPlus->quantumNumbers();
  std::cout << rho->quantumNumbers();
  std::cout << omega->quantumNumbers();*/

  assert(D->consistent());
  D->optimizeSpinAmplitudeSharing();
  assert(D->consistent());

  //std::cout << "rho " << rho->quantumNumbers();
  D->printDecayChain();

  std::cout << "alright! \n";
}
