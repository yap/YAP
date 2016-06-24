#include "logging.h"
#include "BreitWigner.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "ImportanceSampler.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PHSP.h"
#include "Resonance.h"
#include "ZemachFormalism.h"

#include <memory>
#include <random>
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
    D->addChannel({piPlus, piMinus, piPlus});

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices \n";

    // std::cout << "\nHelicity angle symmetrizations with " << M.helicityAngles()->maxSymmetrizationIndex() + 1 << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // get default Dalitz axes
    auto massAxes = M.massAxes();

    // generate point randomly in phase space of model
    std::mt19937 g(0);
    auto P = phsp(M, massAxes, g, 10);

    // create data set with 1 empty data point
    auto data = M.createDataSet(1);

    LOG(INFO) << "BEFORE";
    M.fourMomenta()->printMasses(data[0]);

    LOG(INFO) << "setting squared mass ...";
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        data[0].setFinalStateMomenta(P);
        for (unsigned i = 0; i < 10; ++i)
            data.add(P);
    }

    LOG(INFO) << "AFTER";
    M.fourMomenta()->printMasses(data[0]);

    for (auto p : M.fourMomenta()->finalStateMomenta(data[0]))
        LOG(INFO) << p;

    M.calculate(data);
    auto A_DT = amplitude(M.initialStateParticle()->decayTrees(), data[0]);
    LOG(INFO) << "A_DT = " << A_DT;
    LOG(INFO) << "|A_DT|^2 = " << norm(A_DT);

    yap::ModelIntegral MI(M.initialStateParticle()->decayTrees().at(0));
    yap::ImportanceSampler::calculate(MI, data);

    // for (const auto& kv : MI.diagonals())
    //     LOG(INFO) << to_string(kv.second);
    // for (const auto& kv : MI.offDiagonals())
    //     LOG(INFO) << to_string(kv.second);
    LOG(INFO) << "integral = " << to_string(MI.integral());
    auto ff = fit_fractions(MI);
    for (size_t i = 0; i < ff.size(); ++i)
        LOG(INFO) << "fit fraction DT " << i << " = " << ff[i];
    LOG(INFO) << "sum of fit fractions = " << std::accumulate(ff.begin(), ff.end(), 0.);

    LOG(INFO) << "cached integral components:";
    auto I_cached = cached_integrals(MI);
    for (const auto& row : I_cached)
        LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                     [](const std::string & s, const std::complex<double>& c)
    { return s + "\t" + std::to_string(real(c)) + " + " + std::to_string(imag(c));}).erase(0, 1);

    LOG(INFO) << "integral components:";
    auto I = integrals(MI);
    for (const auto& row : I)
        LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                     [](const std::string & s, const std::complex<double>& c)
    { return s + "\t" + std::to_string(real(c)) + " + " + std::to_string(imag(c));}).erase(0, 1);

    LOG(INFO) << M.initialStateParticle()->printDecayTrees();

    LOG(INFO) << "alright!";
}
