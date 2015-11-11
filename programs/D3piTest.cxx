#include "logging.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "ParticleFactory.h"
#include "Resonance.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <assert.h>
#include <memory>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);
    
    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // final state particles
    auto piPlus = factory.createFinalStateParticle(211, {0, 2});
    auto piMinus = factory.createFinalStateParticle(-211, {1});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // use only L up to 4
    unsigned max2L(2 * 4);

    // rho
    auto rho = std::make_shared<yap::Resonance>(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(rho->massShape()).width()->setValue(0.149);
    rho->addChannels(piPlus, piMinus, max2L);

    // f_2(1270)
    auto f_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_2->massShape()).width()->setValue(0.185);
    f_2->addChannels(piPlus, piMinus, max2L);
    
    // f_0(980)
    auto f_0_980 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_980->massShape()).width()->setValue(0.329);
    f_0_980->addChannels(piPlus, piMinus, max2L);

    // f_0(1370)
    auto f_0_1370 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1370->massShape()).width()->setValue(0.250);
    f_0_1370->addChannels(piPlus, piMinus, max2L);

    // f_0(1500)
    auto f_0_1500 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1500->massShape()).width()->setValue(0.109);
    f_0_1500->addChannels(piPlus, piMinus, max2L);

    // sigma a.k.a. f_0(500)
    auto sigma = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(sigma->massShape()).width()->setValue(0.800);
    sigma->addChannels(piPlus, piMinus, max2L);

    // initial state particle
    auto D = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    
    // Add channels to D
    D->addChannels(rho,      piPlus, max2L);
    D->addChannels(f_2,      piPlus, max2L);
    D->addChannels(f_0_980,  piPlus, max2L);
    D->addChannels(f_0_1370, piPlus, max2L);
    D->addChannels(f_0_1500, piPlus, max2L);
    D->addChannels(sigma,    piPlus, max2L);
    
    // consistency and optimizations
    assert(D->prepare());
    std::cout << "consistent! \n";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";

    std::cout << "\nHelicity angle symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    D->printSpinAmplitudes();
    D->printDataAccessors(false);

    std::cout << "alright! \n";
}
