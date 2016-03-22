#include "logging.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "ZemachFormalism.h"

#include <memory>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    LOG(INFO) << "Start";

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    LOG(INFO) << "Model created";

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    LOG(INFO) << "factory created";

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    LOG(INFO) << "D created";

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    LOG(INFO) << "fsp's created";

    // set final state
    M.setFinalState({piPlus, piMinus, piPlus});

    LOG(INFO) << "final state set";

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

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->maxSymmetrizationIndex() + 1 << " indices \n";

    std::cout << "\nHelicity angle symmetrizations with " << M.helicityAngles()->maxSymmetrizationIndex() + 1 << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

    std::vector<double> m2 = {1, 1};//{0.9, 1.1}; //{0.1, 4};

    // create data set with 1 empty data point
    auto data = M.createDataSet(1);

    LOG(INFO) << "BEFORE";
    M.fourMomenta()->printMasses(data[0]);

    LOG(INFO) << "setting squared mass ...";
    auto P = M.calculateFourMomenta(massAxes, m2);
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        data[0].setFinalStateMomenta(P);
    }

    LOG(INFO) << "AFTER";
    M.fourMomenta()->printMasses(data[0]);

    for (auto p : M.fourMomenta()->finalStateMomenta(data[0]))
        LOG(INFO) << p;

    auto A = M.amplitude(data[0], data);
    LOG(INFO) << "A = " << A;

    LOG(INFO) << "alright!";
}
