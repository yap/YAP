// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "D_K0pi0pi0.h"

#include <BreitWigner.h>
#include <Constants.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <ParticleCombination.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <Resonance.h>

#include <complex>

// ---------------------------------------------------------
std::unique_ptr<yap::Model> D_K0pi0pi0(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
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
    D->addChannel({sigma, Kshort})->freeAmplitudes()[0]->setValue(std::polar(0.67, yap::rad(140.)));
    // D->addChannel({sigma, Kshort})->freeAmplitudes()[0]->setValue(std::polar(0.99, yap::rad(39.)));

    // f_0(980) (as Flatte)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.406, piZero->mass()->value());
    f_0_980_flatte->addChannel(0.406 * 2, Kshort->mass()->value());
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addChannel({piZero, piZero});
    D->addChannel({f_0_980, Kshort})->freeAmplitudes()[0]->setValue(std::polar(1.71, yap::rad(35.2)));
    // D->addChannel({f_0_980, Kshort})->freeAmplitudes()[0]->setValue(std::polar(2.59, yap::rad(44.8)));

    // f_0(1370)
    auto f_0_1370 = F.resonance(30221, radialSize, std::make_shared<yap::BreitWigner>());
    f_0_1370->addChannel({piZero, piZero});
    D->addChannel({f_0_1370, Kshort})->freeAmplitudes()[0]->setValue(std::polar(5.72, yap::rad(340.3)));
    // D->addChannel({f_0_1370, Kshort})->freeAmplitudes()[0]->setValue(std::polar(11.6, yap::rad(15.8)));

    // f_0(1500)
    // auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, std::make_shared<yap::BreitWigner>());
    // f_0_1500->addChannel({piZero, piZero});
    // D->addChannel({f_0_1500, Kshort})->freeAmplitudes()[0]->setValue(Complex_0);
    // D->addChannel({f_0_1500, Kshort})->freeAmplitudes()[0]->setValue(std::polar(20.9, yap::rad(281.4)));

    // f_2(1270)
    auto f_2_1270 = F.resonance(225, radialSize, std::make_shared<yap::BreitWigner>());
    f_2_1270->addChannel({piZero, piZero});
    D->addChannel({f_2_1270, Kshort})->freeAmplitudes()[0]->setValue(std::polar(1.57, yap::rad(282.)));
    // D->addChannel({f_2_1270, Kshort})->freeAmplitudes()[0]->setValue(std::polar(1.57, yap::rad(282.)));

    // K*(892)
    auto Kstar_892 = F.resonance(313, radialSize, std::make_shared<yap::BreitWigner>());
    Kstar_892->addChannel({Kshort, piZero});
    D->addChannel({Kstar_892, piZero})->freeAmplitudes()[0]->setValue(std::polar(1., 0.));

    // K*_2(1430)
    auto Kstar_2_1430 = F.resonance(315, radialSize, std::make_shared<yap::BreitWigner>());
    Kstar_2_1430->addChannel({Kshort, piZero});
    D->addChannel({Kstar_2_1430, piZero})->freeAmplitudes()[0]->setValue(std::polar(0.43, yap::rad(141.)));
    // D->addChannel({Kstar_2_1430, piZero})->freeAmplitudes()[0]->setValue(std::polar(0.85, yap::rad(159.)));

    // K*(1680)
    auto Kstar_1680 = F.resonance(30313, radialSize, std::make_shared<yap::BreitWigner>());
    Kstar_1680->addChannel({Kshort, piZero});
    D->addChannel({Kstar_1680, piZero})->freeAmplitudes()[0]->setValue(std::polar(5.65, yap::rad(55.)));
    // D->addChannel({Kstar_1680, piZero})->freeAmplitudes()[0]->setValue(std::polar(7.07, yap::rad(18.7)));

    return M;
}
