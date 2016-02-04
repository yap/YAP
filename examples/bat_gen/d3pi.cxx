// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "d3pi.h"

#include <BAT/BCMath.h>

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"

#include <complex>
#include <limits>

// ---------------------------------------------------------
d3pi::d3pi(std::string name)
    : BCModel(name)
{
    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(211);
    std::shared_ptr<yap::FinalStateParticle> piMinus = factory.createFinalStateParticle(-211);

    // initial state particle
    D_ = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    D_->setFinalStateParticles({piPlus, piMinus, piPlus});

    // rho
    std::shared_ptr<yap::Resonance> rho = std::make_shared<yap::Resonance>(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(rho->massShape()).width()->setValue(0.149);
    rho->addChannel({piPlus, piMinus});

    // // f_2(1270)
    std::shared_ptr<yap::Resonance> f_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_2->massShape()).width()->setValue(0.185);
    f_2->addChannel({piPlus, piMinus});
    /*
        // f_0(980)
        std::shared_ptr<yap::Resonance> f_0_980 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(f_0_980->massShape()).width()->setValue(0.329);
        f_0_980->addChannel({piPlus, piMinus});

        // f_0(1370)
        std::shared_ptr<yap::Resonance> f_0_1370 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(f_0_1370->massShape()).width()->setValue(0.250);
        f_0_1370->addChannel({piPlus, piMinus});

        // f_0(1500)
        std::shared_ptr<yap::Resonance> f_0_1500 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(f_0_1500->massShape()).width()->setValue(0.109);
        f_0_1500->addChannel({piPlus, piMinus});

        // sigma a.k.a. f_0(500)
        std::shared_ptr<yap::Resonance> sigma = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(sigma->massShape()).width()->setValue(0.800);
        sigma->addChannel({piPlus, piMinus});

        // f_0(500+100i)
        std::shared_ptr<yap::Resonance> f_0_500_100 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), .500, "f_0_500_100", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(f_0_500_100->massShape()).width()->setValue(0.100);
        f_0_500_100->addChannel({piPlus, piMinus});

        // f_0(1500+100i)
        std::shared_ptr<yap::Resonance> f_0_1500_100 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.500, "f_0_1500_100", radialSize, std::make_unique<yap::BreitWigner>());
        static_cast<yap::BreitWigner&>(f_0_1500_100->massShape()).width()->setValue(0.100);
        f_0_1500_100->addChannel({piPlus, piMinus});

    */
    // Add channels to D
    D_->addChannel({rho,      piPlus});
    D_->addChannel({f_2,      piPlus});
    // \todo do not store/prune DataAccessors from InitialStateParticle that are not used
    // D_->addChannel({f_0_980,  piPlus});
    // D_->addChannel({f_0_1370, piPlus});
    // D_->addChannel({f_0_1500, piPlus});
    // D_->addChannel({sigma,    piPlus});

    // D_->addChannel({f_0_500_100,  piPlus});
    // D_->addChannel({f_0_1500_100, piPlus});

    D_->prepare();

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    // unsigned i = 0;
    // freeAmps[i++]->setValue(std::polar(1.,     0.)); // rho
    // freeAmps[i++]->setValue(std::polar(2.1, -123. * TMath::Pi() / 180.)); // f_2
    // freeAmps[i++]->setValue(std::polar(1.4,   12. * TMath::Pi() / 180.)); // f_0_980
    // freeAmps[i++]->setValue(std::polar(1.3,  -21. * TMath::Pi() / 180.)); // f_0_1370
    // freeAmps[i++]->setValue(std::polar(1.1,  -44. * TMath::Pi() / 180.)); // f_0_1500
    // freeAmps[i++]->setValue(std::polar(3.7,   -3. * TMath::Pi() / 180.)); // sigma

    D_->initializeForMonteCarloGeneration(4);
    std::cout << "number of data partitions = " << D_->dataPartitions().size() << std::endl;

    MassAxes_ = D_->getMassAxes({{0, 1}, {1, 2}});

    for (auto& pc : MassAxes_) {
        std::string axis_label = "m2_" + indices_string(*pc).substr(1, 2);
        auto mrange = D_->getMassRange(pc);
        AddParameter(axis_label, pow(mrange[0], 2), pow(mrange[1], 2));
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << pow(mrange[0], 2) << ", " << pow(mrange[1], 2) << "]"
                  << std::endl;
    }

}

// ---------------------------------------------------------
double d3pi::LogLikelihood(const std::vector<double>& parameters)
{
    // unsigned c = GetCurrentChain();
    // return D_->logOfSquaredAmplitude(D_->dataSet()[c], c);
    return 0;
}

// ---------------------------------------------------------
double d3pi::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = D_->calculateFourMomenta(MassAxes_, parameters);

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();
    return 0;

    unsigned c = GetCurrentChain();
    D_->setFinalStateFourMomenta(D_->dataSet()[c], P, c);
    return 0;
}
