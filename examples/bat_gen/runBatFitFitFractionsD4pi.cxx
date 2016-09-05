// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d4pi.h"

#include <DecayTree.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <PHSP.h>
#include <ZemachFormalism.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <random>

const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }

int main()
{
    yap::plainLogs(el::Level::Info);

    // create bat_fit object
    fit_fitFraction m("D4PI_frac_fit", d4pi());

    double D_mass = 1.8648400;

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
    for (auto& dt : decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(0)))
        m.setFitFraction(dt, 43.3e-2,   quad(2.5e-2, 1.9e-2));
    for (auto& dt : decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(1)))
            m.setFitFraction(dt, 2.5e-2,    quad(0.5e-2, 0.4e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), 8.3e-2,    quad(0.7e-2, 0.6e-2));
    /// \todo needed in canonical basis, values are in transversity basis
    amplitude_basis::canonical<double> c(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(std::complex<double>( 1.1e-2),  quad(0.3e-2, 0.3e-2)),
            complex_basis::cartesian<double>(std::complex<double>( 6.4e-2),  quad(0.6e-2, 0.5e-2)),
            complex_basis::cartesian<double>(std::complex<double>(16.88e-2), quad(1.0e-2, 0.8e-2))));

    for (unsigned l = 0; l<3; ++l)
        for (auto& dt : decay_trees(*D, yap::from(D), yap::to(rho), yap::l_equals(l)))
            m.setFitFraction(dt, real(c.amplitudes()[l]), c.covariance()[l][l][0][0]);

    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_0)),   2.4e-2,  quad(2.4e-2, 0.4e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_2)),   4.9e-2,  quad(4.9e-2, 0.5e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(sigma)), 8.2e-2,  quad(8.2e-2, 0.7e-2));

    // set free amplitude parameters of fit
    for (auto& fa : free_amplitudes(*D, yap::from(D), yap::to(a_1)))
        m.fix(fa, 1., 0.);
    //for (auto& fa : free_amplitudes(*D, yap::from(a_1), yap::to(rho), yap::l_equals(0)))
    //    m.fix(fa, 1., 0.);
    for (auto& fa : free_amplitudes(*D, yap::from(a_1), yap::to(rho), yap::l_equals(1)))
        m.fix(fa, 0., 0.);
    for (auto& fa : free_amplitudes(*D, yap::from(a_1), yap::to(rho), yap::l_equals(2)))
        m.setPrior(fa, new BCGaussianPrior(0.241, quad(0.033, 0.024)), new BCGaussianPrior( 82., quad(5.,   4.)));

    m.setPrior(free_amplitude(*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), new BCGaussianPrior(0.493, quad(0.026, 0.021)), new BCGaussianPrior(193., quad(4.,   4.)));

    // polar -> cartesian; transversity -> canonical
    amplitude_basis::canonical<double> can(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.624, rad(357.), {quad(0.023, 0.015), quad(3., 3.)})), // A_longitudinal
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.157, rad(120.), {quad(0.027, 0.020), quad(7., 8.)})), // A_parallel
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.384, rad(163.), {quad(0.020, 0.015), quad(3., 3.)})))); // A_perpendicular

    for (unsigned l = 0; l<3; ++l) {
        // cartesian -> polar
        complex_basis::polar<double> polar(can[l]);

        for (auto& fa : free_amplitudes(*D, yap::from(D), yap::to(rho), yap::l_equals(l)))
            m.setPrior(fa, new BCGaussianPrior(polar.value()[0], polar.covariance()[0][0]), new BCGaussianPrior(polar.value()[1], polar.covariance()[1][1]));
    }


    m.setPrior(free_amplitude(*D, yap::from(D), yap::to(f_0)),   new BCGaussianPrior(0.233, quad(0.019, 0.015)), new BCGaussianPrior(261., quad(7., 3.)));
    m.setPrior(free_amplitude(*D, yap::from(D), yap::to(f_2)),   new BCGaussianPrior(0.338, quad(0.021, 0.016)), new BCGaussianPrior(317., quad(4., 4.)));
    m.setPrior(free_amplitude(*D, yap::from(D), yap::to(sigma)), new BCGaussianPrior(0.432, quad(0.027, 0.022)), new BCGaussianPrior(254., quad(4., 5.)));

    LOG(INFO) << "Parameters summary";
    m.GetParameters().PrintSummary();
    m.GetObservables().PrintSummary();

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()),
                              D_mass, m.axes(), m2r, g,
                              std::numeric_limits<unsigned>::max()));
    LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";
    // partition integration data
    m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 2);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e6);
    m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(2e4 / m.GetNChains()));

    m.WriteMarkovChain("output/" + m.GetSafeName() + "_mcmc.root", "RECREATE");

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    for (size_t i = 0; i < D->decayTrees().at(0).size(); ++i)
        std::cout << i << " = " << to_string(*D->decayTrees().at(0)[i]) << std::endl;

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    m.PrintKnowledgeUpdatePlots("output/" + m.GetSafeName() + "_update.pdf", 2, 2, true);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
