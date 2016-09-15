// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include "../bat_fit.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include <AmplitudeBasis.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <complex>
#include <memory>

using namespace yap;

inline std::unique_ptr<Model> d4pi()
{
    auto F = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);
    
    //
    // resonant particles
    //
    
    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<RelativisticBreitWigner>());
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
    for (auto& f : free_amplitudes(*a_1, to(rho), l_equals(0))) {
        /// \todo has 3 amplitudes for 3 spin projections
        // S wave
        *f = 1;
        f->variableStatus() = VariableStatus::fixed;
    }
    for (auto& f : free_amplitudes(*a_1, to(rho), l_equals(1))) {
        /// \todo has 3 amplitudes for 3 spin projections
        // P wave
        *f = 0.;
        f->variableStatus() = VariableStatus::fixed;
    }
    for (auto& f : free_amplitudes(*a_1, to(rho), l_equals(2))) {
        /// \todo has 3 amplitudes for 3 spin projections
        // D wave
        *f = std::polar(0.241, rad(82.));
    }

    // a_1 -> sigma pi 
    for (auto& f : free_amplitudes(*a_1, to(sigma)))
        *f = std::polar(0.439, rad(193.));
    
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
    auto pipiFlat = F.decayingParticle(F.pdgCode("pi0"), radialSize); // just need any spin0 particle
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
        *freeAmp = std::complex<double>(c[l]);
    }
    
    LOG(INFO) << "D Decay trees:";
    LOG(INFO) << to_string(D->decayTrees());

    return M;
}

inline bat_fit d4pi_fit(std::string name, std::vector<std::vector<unsigned> > pcs = {})
{
    bat_fit m(name, d4pi(), pcs);

    for (auto fa : m.freeAmplitudes()) {
        if (fa->variableStatus() == VariableStatus::fixed)
            m.fix(fa, abs(fa->value()), deg(arg(fa->value())));
        /*else {
            m.setPriors(fa, new ConstantPrior(0, 1.5), new ConstantPrior(-180, 180));
            m.setRealImagRanges(fa, -1.5, 1.5, -1.5, 1.5);
            m.setAbsArgRanges(fa, 0, 1.5, -180, 180);
        }*/
    }

    return m;
}

inline fit_fitFraction d4pi_fit_fitFraction()
{
    // create bat_fit object
    fit_fitFraction m("D4PI_frac_fit", d4pi());

    //double D_mass = 1.8648400;

    m.GetParameter("N_1").Fix(1);

    // find particles
    auto D     = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D0")));
    auto rho   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("rho0")));
    auto sigma = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(500)")));
    auto a_1   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("a_1+")));
    auto f_0   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0")));
    auto f_2   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_2")));

    LOG(INFO) << m.model()->initialStateParticles().at(D).begin()->first;

    // set fit fractions to fit
    /// \todo does not yet work with more than one decayTree
    //m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(0)), 43.3e-2,   quad(2.5e-2, 1.9e-2));
    //m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(1)), 2.5e-2,    quad(0.5e-2, 0.4e-2));
    m.setFitFraction(decay_tree (*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), 8.3e-2,    quad(0.7e-2, 0.6e-2));

    amplitude_basis::canonical<double> c(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(std::complex<double>( 1.1e-2),  quad(0.3e-2, 0.3e-2)),
            complex_basis::cartesian<double>(std::complex<double>( 6.4e-2),  quad(0.6e-2, 0.5e-2)),
            complex_basis::cartesian<double>(std::complex<double>(16.88e-2), quad(1.0e-2, 0.8e-2))));

    /// \todo does not yet work with more than one decayTree
    //for (unsigned l = 0; l<3; ++l)
    //    m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(rho), yap::l_equals(l)), real(c.amplitudes()[l]), c.covariance()[l][l][0][0]);

    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_0)),   2.4e-2,  quad(2.4e-2, 0.4e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_2)),   4.9e-2,  quad(4.9e-2, 0.5e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(sigma)), 8.2e-2,  quad(8.2e-2, 0.7e-2));

    // set free amplitude parameters of fit
    m.fix(free_amplitude(*D, yap::from(D), yap::to(a_1)), 1., 0.);
    m.fix(free_amplitude(*D, yap::from(a_1), yap::to(rho), yap::l_equals(1)), 0., 0.);
    m.setPriors(free_amplitude(*D, yap::from(a_1), yap::to(rho), yap::l_equals(2)), new BCGaussianPrior(0.241, quad(0.033, 0.024)), new BCGaussianPrior( 82., quad(5.,   4.)));

    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), new BCGaussianPrior(0.493, quad(0.026, 0.021)), new BCGaussianPrior(193., quad(4.,   4.)));

    // polar -> cartesian; transversity -> canonical
    amplitude_basis::canonical<double> can(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.624, rad(357.), {quad(0.023, 0.015), quad(3., 3.)})), // A_longitudinal
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.157, rad(120.), {quad(0.027, 0.020), quad(7., 8.)})), // A_parallel
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.384, rad(163.), {quad(0.020, 0.015), quad(3., 3.)})))); // A_perpendicular

    for (unsigned l = 0; l<3; ++l) {
        // cartesian -> polar
        complex_basis::polar<double> polar(can[l]);

        m.setPriors(free_amplitude(*D, yap::from(D), yap::to(rho), yap::l_equals(l)),
                new BCGaussianPrior(polar.value()[0], polar.covariance()[0][0]),
                new BCGaussianPrior(polar.value()[1], polar.covariance()[1][1]));
    }


    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_0)),   new BCGaussianPrior(0.233, quad(0.019, 0.015)), new BCGaussianPrior(261., quad(7., 3.)));
    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_2)),   new BCGaussianPrior(0.338, quad(0.021, 0.016)), new BCGaussianPrior(317., quad(4., 4.)));
    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(sigma)), new BCGaussianPrior(0.432, quad(0.027, 0.022)), new BCGaussianPrior(254., quad(4., 5.)));

    return m;
}

#endif
