#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "HelicitySpinAmplitude.h"
#include "InitialStateParticle.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "SpinUtilities.h"
#include "WignerD.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <assert.h>
#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP


int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);


    /// \todo Figure out clever way to find PDL file
    yap::ParticleFactory factory("evt.pdl");

    // initial state particle
    double radialSize = 1.;
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(421, radialSize);

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(211, {0, 2});
    std::shared_ptr<yap::FinalStateParticle> piMinus = factory.createFinalStateParticle(-211, {1, 3});

    // rho rho
    std::shared_ptr<yap::Resonance> rho = factory.createResonanceBreitWigner(113, radialSize);
    rho->addChannel(piPlus, piMinus, 2 * 1);

    //D->addChannel(rho, rho, 2 * 0); // Clebsch-Gordan coeffs == 0
    //D->addChannel(rho, rho, 2 * 1); // Clebsch-Gordan coeffs == 0
    D->addChannel(rho, rho, 2 * 2);

    // omega omega
    std::shared_ptr<yap::Resonance> omega = factory.createResonanceBreitWigner(223, radialSize);
    omega->addChannel(piPlus, piMinus, 2 * 1);

    //D->addChannel(omega, omega, 2 * 0); // Clebsch-Gordan coeffs == 0
    //D->addChannel(omega, omega, 2 * 1); // Clebsch-Gordan coeffs == 0
    D->addChannel(omega, omega, 2 * 2);

    // rho omega
    //D->addChannel(rho, omega, 2 * 0); // Clebsch-Gordan coeffs == 0
    //D->addChannel(rho, omega, 2 * 1); // Clebsch-Gordan coeffs == 0
    D->addChannel(rho, omega, 2 * 2);


    // a_1 channels
    std::shared_ptr<yap::Resonance> sigma = factory.createResonanceBreitWigner(9000221, radialSize);
    sigma->addChannel(piPlus, piMinus, 2 * 0);

    std::shared_ptr<yap::Resonance> a_1 = factory.createResonanceBreitWigner(20213, radialSize);
    a_1->addChannel(sigma, piPlus, 2 * 1);

    a_1->addChannel(rho, piPlus, 2 * 0); // S-wave
    a_1->addChannel(rho, piPlus, 2 * 1); // not in Focus model
    a_1->addChannel(rho, piPlus, 2 * 2); // D-wave

    D->addChannel(a_1, piMinus, 2 * 1);


    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.createResonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);




    // consistency and optimizations
    assert(D->prepare());
    std::cout << "consistent! \n";

    // print stuff
    yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\nD symmetrizations: \n";
    for (std::shared_ptr<yap::ParticleCombination>& pc : D->particleCombinations())
        std::cout << std::string(*pc) << "\n";
    std::cout << "\n";

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex()+1 << " indices: \n";
    for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex()+1 << " indices: \n";
    for (auto& pc : D->helicityAngles().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->helicityAngles().symmetrizationIndex(pc) << "\n";

    D->printDecayChain();
    std::cout << "\n";

    D->printSpinAmplitudes();
    D->printDataAccessors();



    // test helicity angles
    TLorentzVector P(0.0, 0.0, 0.0, D->mass());
    Double_t masses[4] = { piPlus->mass(), piMinus->mass(), piPlus->mass(), piMinus->mass() };

    TGenPhaseSpace event;
    event.SetDecay(P, 4, masses);
    event.Generate();

    std::vector<TLorentzVector> momenta;
    for (unsigned i = 0; i < 4; ++i)
        momenta.push_back(*event.GetDecay(i));

    D->addDataPoint(yap::DataPoint(momenta));



    std::cout << "alright! \n";
}
