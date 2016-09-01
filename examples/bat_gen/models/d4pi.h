// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include <AmplitudeBasis.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

using namespace yap;

inline std::unique_ptr<Model> d4pi(std::unique_ptr<yap::Model> M)
{
    auto F = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);
    
    //
    // resonant particles
    //
    
    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<BreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // omega
    //auto omega = F.resonance(F.pdgCode("omega"), radialSize, std::make_shared<BreitWigner>());
    //omega->addChannel({piPlus, piMinus});
    
    // sigma / f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<BreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = F.resonance(F.pdgCode("a_1+"), radialSize, std::make_shared<BreitWigner>());
    a_1->addChannel(rho,   piPlus);
    a_1->addChannel(sigma, piPlus);

    // a_1 -> rho pi
    *free_amplitude(*a_1, to(rho), l_equals(0)) = 1; // S wave
    free_amplitude(*a_1, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
    *free_amplitude(*a_1, to(rho), l_equals(1)) = 0.; // P wave
    *free_amplitude(*a_1, to(rho), l_equals(2)) = std::polar(0.241, rad(82.)); // D wave

    // a_1 -> sigma pi 
    *free_amplitude(*a_1, to(sigma)) = std::polar(0.439, rad(193.));
    
    // f_0(980) (as Flatte)
    auto piZero = F.fsp(111);
    auto Kshort = F.fsp(310);
    auto f_0_980_flatte = std::make_shared<Flatte>();
    f_0_980_flatte->add(FlatteChannel(0.20, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.50, *F.fsp(321), *F.fsp(-321))); // K+K-
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
       
    // f_2(1270)
    auto f_2 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<BreitWigner>());
    f_2->addChannel(piPlus, piMinus); 
    
    // pi+ pi- flat
    auto pipiFlat = F.decayingParticle(F.pdgCode("f_0"), radialSize); // just need any spin0 particle
    pipiFlat->addChannel(piPlus, piMinus);   
    
    //
    // D0 channels
    //

    D->addChannel(rho, rho);
    D->addChannel(a_1, piMinus);
    D->addChannel(f_0_980, piPlus, piMinus);
    D->addChannel(f_2, pipiFlat);
    D->addChannel(sigma, piPlus, piMinus);
    
    M->addInitialStateParticle(D);

    // R pi pi
    *free_amplitude(*M, to(f_0_980, piPlus, piMinus)) = std::polar(0.233, rad(261.));
    *free_amplitude(*M, to(f_2,     pipiFlat       )) = std::polar(0.338, rad(317.));
    *free_amplitude(*M, to(sigma,   piPlus, piMinus)) = std::polar(0.432, rad(254.));
    
    // rho rho
    // transform into angular momentum basis
    amplitude_basis::canonical<double> c(amplitude_basis::transversity<double>(
                                             std::polar(0.624, rad(357.)),    // A_longitudinal
                                             std::polar(0.157, rad(120.)),    // A_parallel
                                             std::polar(0.384, rad(163.)) )); // A_perpendicular
    
    for (unsigned l = 0; l < 3; ++l) {
        auto freeAmp = free_amplitude(*M, to(rho, rho), l_equals(l));
        LOG(INFO) << to_string(*freeAmp);
        *freeAmp = c[l];
    }
    
    return M;
}

#endif
