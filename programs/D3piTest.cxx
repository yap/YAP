#include "logging.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"

#include <memory>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]
    // use only L up to 4
    unsigned maxL(2);

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

    // initial state particle
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(211);
    std::shared_ptr<yap::FinalStateParticle> piMinus = factory.createFinalStateParticle(-211);

    // set final state
    D->setFinalStateParticles({piPlus, piMinus, piPlus});

    // rho
    std::shared_ptr<yap::Resonance> rho = std::make_shared<yap::Resonance>(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(rho->massShape()).width()->setValue(0.149);
    rho->addChannels(piPlus, piMinus, maxL);

    // f_2(1270)
    std::shared_ptr<yap::Resonance> f_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_2->massShape()).width()->setValue(0.185);
    f_2->addChannels(piPlus, piMinus, maxL);

    // f_0(980)
    std::shared_ptr<yap::Resonance> f_0_980 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_980->massShape()).width()->setValue(0.329);
    f_0_980->addChannels(piPlus, piMinus, maxL);

    // f_0(1370)
    std::shared_ptr<yap::Resonance> f_0_1370 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1370->massShape()).width()->setValue(0.250);
    f_0_1370->addChannels(piPlus, piMinus, maxL);

    // f_0(1500)
    std::shared_ptr<yap::Resonance> f_0_1500 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1500->massShape()).width()->setValue(0.109);
    f_0_1500->addChannels(piPlus, piMinus, maxL);

    // sigma a.k.a. f_0(500)
    std::shared_ptr<yap::Resonance> sigma = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(sigma->massShape()).width()->setValue(0.800);
    sigma->addChannels(piPlus, piMinus, maxL);

    // Add channels to D
    D->addChannels(rho,      piPlus, maxL);
    D->addChannels(f_2,      piPlus, maxL);
    D->addChannels(f_0_980,  piPlus, maxL);
    D->addChannels(f_0_1370, piPlus, maxL);
    D->addChannels(f_0_1500, piPlus, maxL);
    D->addChannels(sigma,    piPlus, maxL);

    // consistency and optimizations
    D->prepare();
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

    // initialize for 5 streams
    D->initializeForMonteCarloGeneration(5);

    // choose Dalitz coordinates m^2_12 and m^2_23
    yap::ParticleCombinationVector DalitzAxes = D->fourMomenta().getDalitzAxes({{0, 1}, {1, 2}});

    std::vector<double> m2(DalitzAxes.size(), 1);

    DEBUG("BEFORE");
    D->fourMomenta().printMasses(D->dataSet()[0]);

    D->fourMomenta().setSquaredMasses(D->dataSet()[0], DalitzAxes, m2);

    DEBUG("AFTER");
    D->fourMomenta().printMasses(D->dataSet()[0]);

    std::cout << "alright! \n";
}
