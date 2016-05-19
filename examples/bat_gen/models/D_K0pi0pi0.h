// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKSPI0PI0__H
#define __BAT__DKSPI0PI0__H


#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <Resonance.h>
#include <Parameter.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>


inline std::unique_ptr<yap::Model> D_K0pi0pi0(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{
    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piZero = F.fsp(F.pdgCode("pi0"));
    auto Kshort = F.fsp(F.pdgCode("K_S0"));

    auto M = std::make_unique<yap::Model>(std::move(SAC));
    M->setFinalState({Kshort, piZero, piZero});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);

    // f_0(500), aka "sigma" (as PoleMass)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<yap::PoleMass>(std::complex<double>(0.470, 0.220)));
    sigma->addChannel({piZero, piZero});

    // f_0(980) (as Flatte)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.406, piZero->mass()->value());
    f_0_980_flatte->addChannel(0.406 * 2, Kshort->mass()->value());
    auto f_0_980 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.965, "f_0_980", radialSize, f_0_980_flatte);
    f_0_980->addChannel({piZero, piZero});

    // f_0(1370)
    auto f_0_1370 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 1.350, "f_0_1370", radialSize, std::make_shared<yap::BreitWigner>(0.265));
    f_0_1370->addChannel({piZero, piZero});

    // f_0(1500)
    auto f_0_1500 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 1.505, "f_0_1500", radialSize, std::make_shared<yap::BreitWigner>(0.109));
    f_0_1500->addChannel({piZero, piZero});

    // f_2(1270)
    auto f_2_1270 = yap::Resonance::create(yap::QuantumNumbers(2 * 2, 0), 1.2751, "f_2_1270", radialSize, std::make_shared<yap::BreitWigner>(0.185));
    f_2_1270->addChannel({piZero, piZero});

    // K*(892)
    auto Kstar_892 = yap::Resonance::create(yap::QuantumNumbers(1 * 2, 0), 0.896, "Kstar_892", radialSize, std::make_shared<yap::BreitWigner>(0.0503));
    Kstar_892->addChannel({Kshort, piZero});

    // K*_2(1430)
    auto Kstar_2_1430 = yap::Resonance::create(yap::QuantumNumbers(2 * 2, 0), 1.4324, "Kstar_2_1430", radialSize, std::make_shared<yap::BreitWigner>(0.109));
    Kstar_2_1430->addChannel({Kshort, piZero});

    // K*(1680)
    auto Kstar_1680 = yap::Resonance::create(yap::QuantumNumbers(1 * 2, 0), 1.717, "Kstar_1680", radialSize, std::make_shared<yap::BreitWigner>(0.322));
    Kstar_1680->addChannel({Kshort, piZero});

    // change switch argument to choose between different CLEO models
    switch (1) {

        case 1:

            D->addChannel({sigma, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(0.67, yap::rad(140.)));
            D->addChannel({f_0_980, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(1.71, yap::rad(35.2)));
            D->addChannel({f_0_1370, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(5.72, yap::rad(340.3)));
            D->addChannel({f_0_1500, Kshort})->freeAmplitudes().begin()->get()->setValue(yap::Complex_0);
            D->addChannel({f_2_1270, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(1.57, yap::rad(282.)));
            D->addChannel({Kstar_892, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(1., 0.));
            D->addChannel({Kstar_2_1430, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(0.43, yap::rad(141.)));
            D->addChannel({Kstar_1680, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(5.65, yap::rad(55.)));

            break;

        case 2:

            D->addChannel({sigma, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(0.91, yap::rad(119.)));
            D->addChannel({f_0_980, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(2.13, yap::rad(65.)));
            D->addChannel({f_0_1370, Kshort})->freeAmplitudes().begin()->get()->setValue(yap::Complex_0);
            D->addChannel({f_0_1500, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(11.7, yap::rad(16.)));
            D->addChannel({f_2_1270, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(4.16, yap::rad(2.2)));
            D->addChannel({Kstar_892, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(1., 0.));
            D->addChannel({Kstar_2_1430, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(0.98, yap::rad(191.)));
            D->addChannel({Kstar_1680, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(7.6, yap::rad(45.)));

            break;

        case 3:

            D->addChannel({sigma, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(0.99, yap::rad(39.)));
            D->addChannel({f_0_980, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(2.59, yap::rad(44.8)));
            D->addChannel({f_0_1370, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(11.6, yap::rad(15.8)));
            D->addChannel({f_0_1500, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(20.9, yap::rad(281.4)));
            D->addChannel({f_2_1270, Kshort})->freeAmplitudes().begin()->get()->setValue(std::polar(2.98, yap::rad(340.9)));
            D->addChannel({Kstar_892, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(1., 0.));
            D->addChannel({Kstar_2_1430, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(0.85, yap::rad(159.)));
            D->addChannel({Kstar_1680, piZero})->freeAmplitudes().begin()->get()->setValue(std::polar(7.07, yap::rad(18.7)));

            break;

        default:
            throw std::runtime_error("CLEO Model must be 1, 2, or 3");
    }

    return M;
}

#endif
