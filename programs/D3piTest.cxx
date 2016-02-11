#include "logging.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "HelicitySpinAmplitude.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "ZemachSpinAmplitude.h"

#include <memory>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    LOG(DEBUG) << "1";

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    LOG(DEBUG) << "2";

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    LOG(DEBUG) << "3";

    // set final state
    M.setFinalState({piPlus, piMinus, piPlus});

    LOG(DEBUG) << "4";

    // rho
    auto rho = yap::Resonance::create(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(rho->massShape())->width()->setValue(0.149);
    rho->addChannel({piPlus, piMinus});

    // f_2(1270)
    auto f_2 = yap::Resonance::create(factory.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_2->massShape())->width()->setValue(0.185);
    f_2->addChannel({piPlus, piMinus});

    // f_0(980)
    auto f_0_980 = yap::Resonance::create(factory.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_980->massShape())->width()->setValue(0.329);
    f_0_980->addChannel({piPlus, piMinus});

    // f_0(1370)
    auto f_0_1370 = yap::Resonance::create(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->width()->setValue(0.250);
    f_0_1370->addChannel({piPlus, piMinus});

    // f_0(1500)
    auto f_0_1500 = yap::Resonance::create(factory.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_1500->massShape())->width()->setValue(0.109);
    f_0_1500->addChannel({piPlus, piMinus});

    // sigma a.k.a. f_0(500)
    auto sigma = yap::Resonance::create(factory.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(sigma->massShape())->width()->setValue(0.800);
    sigma->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});
    D->addChannel({f_2,      piPlus});
    D->addChannel({f_0_980,  piPlus});
    D->addChannel({f_0_1370, piPlus});
    D->addChannel({f_0_1500, piPlus});
    D->addChannel({sigma,    piPlus});

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";

    std::cout << "\nHelicity angle symmetrizations with " << M.helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // initialize for 5 streams
    M.initializeForMonteCarloGeneration(5);

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.getMassAxes({{0, 1}, {1, 2}});

    std::vector<double> m2 = {1, 1};//{0.9, 1.1}; //{0.1, 4};

    LOG(INFO) << "BEFORE";
    M.fourMomenta().printMasses(M.dataSet()[0]);

    LOG(INFO) << "setting squared mass ...";
    auto P = M.calculateFourMomenta(massAxes, m2);
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        M.setFinalStateMomenta(M.dataSet()[0], P);
    }

    LOG(INFO) << "AFTER";
    M.fourMomenta().printMasses(M.dataSet()[0]);

    for (auto p : M.fourMomenta().finalStateMomenta(M.dataSet()[0]))
        LOG(INFO) << p;

    M.resetCalculationStatuses(0);
    auto A = M.amplitude(M.dataSet()[0], 0);
    LOG(INFO) << "A = " << A;

    LOG(INFO) << "alright!";
}
