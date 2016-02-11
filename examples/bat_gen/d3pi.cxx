// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "d3pi.h"

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "ParticleCombination.h"
#include "QuantumNumbers.h"
#include "Resonance.h"

#include <complex>

// ---------------------------------------------------------
std::unique_ptr<yap::Model> d3pi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{
    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::move(SAC));
    M->setFinalState({piPlus, piMinus, piPlus});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // rho
    // auto rho = std::make_shared<yap::Resonance>(yap::QuantumNumbers(2, 0), 0.775, "rho", radialSize, std::make_shared<yap::BreitWigner>());
    auto rho = std::make_shared<yap::Resonance>(yap::QuantumNumbers(2, 0), 1.2, "rho", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(rho->massShape())->width()->setValue(0.1);
    rho->addChannel({piPlus, piMinus});
    D->addChannel({rho, piPlus});

    /*
    // // f_2(1270)
    auto f_2 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_2->massShape())->width()->setValue(0.185);
    f_2->addChannel({piPlus, piMinus});

    // f_0(980)
    auto f_0_980 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_980->massShape())->width()->setValue(0.329);
    f_0_980->addChannel({piPlus, piMinus});

    // f_0(1370)
    auto f_0_1370 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->width()->setValue(0.250);
    f_0_1370->addChannel({piPlus, piMinus});

    // f_0(1500)
    auto f_0_1500 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_1500->massShape())->width()->setValue(0.109);
    f_0_1500->addChannel({piPlus, piMinus});

    // sigma a.k.a. f_0(500)
    auto sigma = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(sigma->massShape())->width()->setValue(0.800);
    sigma->addChannel({piPlus, piMinus});

    // f_0(500+100i)
    auto f_0_500_100 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), .500, "f_0_500_100", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_500_100->massShape())->width()->setValue(0.100);
    f_0_500_100->addChannel({piPlus, piMinus});

    // f_0(1500+100i)
    auto f_0_1500_100 = std::make_shared<yap::Resonance>(F.quantumNumbers("f_0"), 1.500, "f_0_1500_100", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(f_0_1500_100->massShape())->width()->setValue(0.100);
    f_0_1500_100->addChannel({piPlus, piMinus});

    */
    // Add channels to D
    // D->addChannel({f_2,      piPlus});
    // \todo do not store/prune DataAccessors from InitialStateParticle that are not used
    // D->addChannel({f_0_980,  piPlus});
    // D->addChannel({f_0_1370, piPlus});
    // D->addChannel({f_0_1500, piPlus});
    // D->addChannel({sigma,    piPlus});

    // D->addChannel({f_0_500_100,  piPlus});
    // D->addChannel({f_0_1500_100, piPlus});

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = M->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    // unsigned i = 0;
    // freeAmps[i++]->setValue(std::polar(1.,     0.)); // rho
    // freeAmps[i++]->setValue(std::polar(2.1, -123. * TMath::Pi() / 180.)); // f_2
    // freeAmps[i++]->setValue(std::polar(1.4,   12. * TMath::Pi() / 180.)); // f_0_980
    // freeAmps[i++]->setValue(std::polar(1.3,  -21. * TMath::Pi() / 180.)); // f_0_1370
    // freeAmps[i++]->setValue(std::polar(1.1,  -44. * TMath::Pi() / 180.)); // f_0_1500
    // freeAmps[i++]->setValue(std::polar(3.7,   -3. * TMath::Pi() / 180.)); // sigma

    return M;
}
