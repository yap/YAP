// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKSPI0PI0__H
#define __BAT__DKSPI0PI0__H

#include "find_pdl_file.h"

#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <MathUtilities.h>
#include <Model.h>
#include <PDL.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>


inline std::unique_ptr<yap::Model> D_K0pi0pi0(std::unique_ptr<yap::Model> M)
{
    auto T = yap::read_pdl_file(find_pdl_file());

    // final state particles
    auto piZero = FinalStateParticle::create(T["pi0"]);
    auto Kshort = FinalStateParticle::create(T["K_S0"]);

    M->setFinalState(Kshort, piZero, piZero);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D0"], radialSize);

    // f_0(500), aka "sigma" (as PoleMass)
    auto sigma = DecayingParticle::create(T["f_0(500)"], radialSize, std::make_shared<yap::PoleMass>(std::complex<double>(0.470, -0.220)));
    sigma->addStrongDecay(piZero, piZero);
    D->addWeakDecay(sigma, Kshort);
    
    // f_0(980) (as Flatte)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>(0.965);
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406, *piZero, *piZero));
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406 * 4, *Kshort, *Kshort));
    auto f_0_980 = yap::DecayingParticle::create("f_0_980", T["f_0"].quantumNumbers(), radialSize, f_0_980_flatte);
    f_0_980->addStrongDecay(piZero, piZero);
    D->addWeakDecay(f_0_980, Kshort);

    // f_0(1370)
    auto f_0_1370 = yap::DecayingParticle::create("f_0_1370", T["f_0"].quantumNumbers(), radialSize, std::make_shared<yap::BreitWigner>(1.35, 0.265));
    f_0_1370->addStrongDecay(piZero, piZero);
    D->addWeakDecay(f_0_1370, Kshort);

    // f_0(1500)
    auto f_0_1500 = yap::DecayingParticle::create("f_0_1500", T["f_0"].quantumNumbers(), radialSize, std::make_shared<yap::BreitWigner>(1.505, 0.109));
    f_0_1500->addStrongDecay(piZero, piZero);
    D->addWeakDecay(f_0_1500, Kshort);

    // f_2(1270)
    auto f_2 = yap::DecayingParticle::create(T["f_2"], radialSize, std::make_shared<yap::BreitWigner>(1.2751, 0.185));
    f_2->addStrongDecay(piZero, piZero);
    D->addWeakDecay(f_2, Kshort);

    // K*(892)
    auto Kstar_892 = yap::DecayingParticle::create("Kstar_892", T["K*0"].quantumNumbers(), radialSize, std::make_shared<yap::BreitWigner>(0.896, 0.0503));
    Kstar_892->addStrongDecay(Kshort, piZero);
    D->addWeakDecay(Kstar_892, piZero);

    // K*_2(1430)
    auto Kstar_2_1430 = yap::DecayingParticle::create("Kstar_2_1430", T["K_2*0"].quantumNumbers(), radialSize, std::make_shared<yap::BreitWigner>(1.4324, .109));
    Kstar_2_1430->addStrongDecay(Kshort, piZero);
    D->addWeakDecay(Kstar_2_1430, piZero);

    // K*(1680)
    auto Kstar_1680 = yap::DecayingParticle::create("Kstar_1680", T["phi(1680)"].quantumNumbers(), radialSize, std::make_shared<yap::BreitWigner>(1.717, 0.322));
    Kstar_1680->addStrongDecay(Kshort, piZero);
    D->addWeakDecay(Kstar_1680, piZero);

    // change switch argument to choose between different CLEO models
    switch (2) {

        case 1:
            
            *free_amplitude(*M, to(sigma))        = std::polar(0.67, yap::rad(140.));
            *free_amplitude(*M, to(f_0_980))      = std::polar(1.71, yap::rad(35.2));
            *free_amplitude(*M, to(f_0_1370))     = std::polar(5.72, yap::rad(340.3));
            *free_amplitude(*M, to(f_0_1500))     = 0.;
            *free_amplitude(*M, to(f_2))     = std::polar(1.57, yap::rad(282.));
            *free_amplitude(*M, to(Kstar_892))    = std::polar(1., 0.);
            *free_amplitude(*M, to(Kstar_2_1430)) = std::polar(0.43, yap::rad(141.));
            *free_amplitude(*M, to(Kstar_1680))   = std::polar(5.65, yap::rad(55.));

            break;

        case 2:

            *free_amplitude(*M, to(sigma))         = std::polar(0.91, yap::rad(119.));
            *free_amplitude(*M, to(f_0_980))       = std::polar(2.13, yap::rad(65.));
            *free_amplitude(*M, to(f_0_1370))      = 0.;
            *free_amplitude(*M, to(f_0_1500))      = std::polar(11.7, yap::rad(16.));
            *free_amplitude(*M, to(f_2))      = std::polar(4.16, yap::rad(2.2));
            *free_amplitude(*M, to(Kstar_892))     = std::polar(1., 0.);
            *free_amplitude(*M, to(Kstar_2_1430))  = std::polar(0.98, yap::rad(191.));
            *free_amplitude(*M, to(Kstar_1680))    = std::polar(7.6, yap::rad(45.));

            break;

        case 3:

            *free_amplitude(*M, to(sigma))        = std::polar(0.99, yap::rad(39.));
            *free_amplitude(*M, to(f_0_980))      = std::polar(2.59, yap::rad(44.8));
            *free_amplitude(*M, to(f_0_1370))     = std::polar(11.6, yap::rad(15.8));
            *free_amplitude(*M, to(f_0_1500))     = std::polar(20.9, yap::rad(281.4));
            *free_amplitude(*M, to(f_2))     = std::polar(2.98, yap::rad(340.9));
            *free_amplitude(*M, to(Kstar_892))    = std::polar(1., 0.);
            *free_amplitude(*M, to(Kstar_2_1430)) = std::polar(0.85, yap::rad(159.));
            *free_amplitude(*M, to(Kstar_1680))   = std::polar(7.07, yap::rad(18.7));

            break;

        default:
            throw std::runtime_error("CLEO Model must be 1, 2, or 3");
    }

    return M;
}

#endif
