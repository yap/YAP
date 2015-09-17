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

    unsigned max2L(2*4);

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
    rho->addChannels(piPlus, piMinus, max2L);

    D->addChannels(rho, rho, max2L);

    // omega omega
    std::shared_ptr<yap::Resonance> omega = factory.createResonanceBreitWigner(223, radialSize);
    omega->addChannels(piPlus, piMinus, max2L);

    D->addChannels(omega, omega, max2L);

    // rho omega
    D->addChannels(rho, omega, max2L);


    // a_1 channels
    std::shared_ptr<yap::Resonance> sigma = factory.createResonanceBreitWigner(9000221, radialSize);
    sigma->addChannels(piPlus, piMinus, max2L);

    std::shared_ptr<yap::Resonance> a_1 = factory.createResonanceBreitWigner(20213, radialSize);
    a_1->addChannels(sigma, piPlus, max2L);

    a_1->addChannels(rho, piPlus, max2L);

    D->addChannels(a_1, piMinus, max2L);


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

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices: \n";
    for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices: \n";
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

    D->logLikelihood();





    std::cout << "alright! \n";
}
