#include "logging.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "MassAxes.h"
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

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

    // create final state particles
    std::shared_ptr<yap::FinalStateParticle> kPlus  = factory.createFinalStateParticle(+321);
    std::shared_ptr<yap::FinalStateParticle> kMinus = factory.createFinalStateParticle(-321);
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(+211);

    // create initial state particle and set final state
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    D->setFinalStateParticles({kPlus, kMinus, piPlus});

    // create a phi
    auto phi = std::make_shared<yap::Resonance>(factory.quantumNumbers("phi"), 1019.461e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(phi->massShape()).width()->setValue(4.266e-3);
    phi->addChannel({kPlus, kMinus});

    // Add channels to D
    D->addChannel({phi, piPlus});

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

    std::cout << D->spinAmplitudeCache() << std::endl;
    D->printDataAccessors(false);

    // initialize for 5 streams
    D->initializeForMonteCarloGeneration(5);

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = D->getMassAxes({{0, 1}, {1, 2}});

    std::vector<double> m2(massAxes.size(), 1);

    DEBUG("BEFORE");
    D->fourMomenta().printMasses(D->dataSet()[0]);

    LOG(INFO) << "setting squared mass ...";
    auto P = D->calculateFourMomenta(massAxes, m2);
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        D->setFinalStateFourMomenta(D->dataSet()[0], P);
    }        

    DEBUG("AFTER");
    D->fourMomenta().printMasses(D->dataSet()[0]);

    std::cout << "alright! \n";
}
