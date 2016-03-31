// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include <BreitWigner.h>
#include <Constants.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <QuantumNumbers.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

inline std::unique_ptr<yap::Model> d3pi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{
    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::move(SAC));
    M->setFinalState({piPlus, piMinus, piPlus});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // f_0(500), aka "sigma"
    auto sigma = yap::Resonance::create(F.quantumNumbers("f_0(500)"), 0.800, "sigma", radialSize, std::make_shared<yap::BreitWigner>(0.800));
    sigma->addChannel({piPlus, piMinus});
    D->addChannel({sigma, piPlus})->freeAmplitudes()[0]->setValue(std::polar(3.7, yap::rad(-3.)));

    // f_0(1370)
    auto f_0_1370 = F.resonance(F.pdgCode("f'_0"), radialSize, std::make_shared<yap::BreitWigner>());
    f_0_1370->addChannel({piPlus, piMinus});
    D->addChannel({f_0_1370, piPlus})->freeAmplitudes()[0]->setValue(std::polar(1.3, yap::rad(-21.)));

    // f_0(1500)
    auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, std::make_shared<yap::BreitWigner>());
    f_0_1500->addChannel({piPlus, piMinus});
    D->addChannel({f_0_1500, piPlus})->freeAmplitudes()[0]->setValue(std::polar(1.1, yap::rad(-44.)));

    // rho^0(770)
    auto rho0 = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<yap::BreitWigner>());
    rho0->addChannel({piPlus, piMinus});
    D->addChannel({rho0, piPlus})->freeAmplitudes()[0]->setValue(yap::Complex_1);

    // f_2(1270)
    auto f_2_1270 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<yap::BreitWigner>());
    f_2_1270->addChannel({piPlus, piMinus});
    D->addChannel({f_2_1270, piPlus})->freeAmplitudes()[0]->setValue(std::polar(2.1, yap::rad(-123.)));

    // // f_0(980) (as Breit-Wigner)
    // auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, std::make_shared<yap::BreitWigner>());
    // f_0_980->addChannel({piPlus, piMinus});
    // D->addChannel({f_0_980, piPlus})->freeAmplitudes()[0]->setValue(std::polar(0.75, yap::rad(12.)));

    // auto pp0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 1.1, "pp0", radialSize, std::make_shared<yap::BreitWigner>(1.1, 0.05));
    // pp0->addChannel({piPlus, piMinus});
    // D->addChannel({pp0, piPlus})->freeAmplitudes()[0]->setValue(yap::Complex_1);

    // auto pp1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.35, "pp1", radialSize, std::make_shared<yap::BreitWigner>(1.35, 0.05));
    // pp1->addChannel({piPlus, piMinus});
    // D->addChannel({pp1, piPlus})->freeAmplitudes()[0]->setValue(2. * yap::Complex_1);

    // auto pp2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.6, "pp2", radialSize, std::make_shared<yap::BreitWigner>(1.6, 0.05));
    // pp2->addChannel({piPlus, piMinus});
    // D->addChannel({pp2, piPlus})->freeAmplitudes()[0]->setValue(30. * yap::Complex_1);

    return M;
}

#endif
