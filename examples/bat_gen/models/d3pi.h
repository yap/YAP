// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <make_unique.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

inline std::unique_ptr<yap::Model> d3pi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::move(SAC));
    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // rho
    auto rho = F.resonance(113, radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    rho->addChannel(piPlus, piMinus);
    D->addChannel(rho, piPlus);

    // f_2(1270)
    auto f_2 = F.resonance(225, radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    f_2->addChannel(piPlus, piMinus);
    D->addChannel(f_2, piPlus);

    // f_0(980)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.406, piPlus->mass()->value());
    f_0_980_flatte->addChannel(0.406 * 2, F.particleTableEntry("K+").Mass);
    auto f_0_980 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.965, "f_0_980", radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
    D->addChannel(f_0_980, piPlus);

    // f_0(1370)
    auto f_0_1370 = yap::Resonance::create(F.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::RelativisticBreitWigner>(0.265));
    f_0_1370->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, std::make_unique<yap::RelativisticBreitWigner>());
    f_0_1500->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_unique<yap::PoleMass>(std::complex<double>(0.470, -0.220)));
    sigma->addChannel(piPlus, piMinus);
    D->addChannel(sigma, piPlus);

    // Add channels to D
    *D->freeAmplitudes(rho,      piPlus)[0] = std::polar(1., 0.);
    *D->freeAmplitudes(f_0_980,  piPlus)[0] = std::polar(1.4, yap::rad(12.));
    *D->freeAmplitudes(f_2,      piPlus)[0] = std::polar(2.1, yap::rad(-123.));
    *D->freeAmplitudes(f_0_1370, piPlus)[0] = std::polar(1.3, yap::rad(-21.));
    *D->freeAmplitudes(f_0_1500, piPlus)[0] = std::polar(1.1, yap::rad(-44.));
    *D->freeAmplitudes(sigma,    piPlus)[0] = std::polar(3.7, yap::rad(-3.));

    M->addInitialStateParticle(D);

    return M;
}

#endif
