// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d3pi.h"

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

    //
    // CREATE yap::Model
    //
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::ZemachFormalism>());
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
    // f_0_980_flatte->addChannel(0.406, piPlus->mass()->value());
    f_0_980_flatte->addChannel(0.329, piPlus->mass()->value());
    f_0_980_flatte->addChannel(2 * f_0_980_flatte->channels().back().Coupling->value(), F.particleTableEntry("K+").Mass);
    // auto f_0_980 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.965, "f_0_980", radialSize, f_0_980_flatte);
    auto f_0_980 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.953, "f_0_980", radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
    D->addChannel(f_0_980, piPlus);

    // f_0(1370)
    // auto f_0_1370 = yap::Resonance::create(F.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::RelativisticBreitWigner>(0.265));
    auto f_0_1370 = yap::Resonance::create(F.quantumNumbers("f_0"), 1.259, "f_0_1370", radialSize, std::make_unique<yap::RelativisticBreitWigner>(0.298));
    f_0_1370->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, std::make_unique<yap::RelativisticBreitWigner>());
    f_0_1500->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    // auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_unique<yap::PoleMass>(std::complex<double>(0.470, -0.220)));
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_unique<yap::PoleMass>(std::complex<double>(0.466, -0.223)));
    sigma->addChannel(piPlus, piMinus);
    D->addChannel(sigma, piPlus);

    M->addInitialStateParticle(D);

    // create bat_fit object
    fit_fitFraction m("D3PI_frac_fit", std::move(M));

    m.GetParameters()[1].Fix(1);

    // set fit fractions to fit
    m.setFitFraction(D->decayTrees(rho,      piPlus)[0], 20e-2,   2.3e-2, 0.9e-2);
    m.setFitFraction(D->decayTrees(f_2,      piPlus)[0], 18.2e-2, 2.6e-2, 0.7e-2);
    m.setFitFraction(D->decayTrees(f_0_980,  piPlus)[0], 4.1e-2,  0.9e-2, 0.3e-2);
    m.setFitFraction(D->decayTrees(f_0_1370, piPlus)[0], 2.6e-2,  1.8e-2, 0.6e-2);
    m.setFitFraction(D->decayTrees(f_0_1500, piPlus)[0], 3.4e-2,  1.0e-2, 0.8e-2);
    m.setFitFraction(D->decayTrees(sigma,    piPlus)[0], 41.8e-2, 1.4e-2, 2.5e-2);

    // set free amplitude parameters of fit
    m.fixFreeAmplitude(m.isp()->freeAmplitudes(rho,      piPlus)[0], 1., 0.);
    m.setFreeAmplitude(m.isp()->freeAmplitudes(f_2,      piPlus)[0], 2.1, quad(0.2, 0.1), -123., quad(6., 3.));
    m.setFreeAmplitude(m.isp()->freeAmplitudes(f_0_980,  piPlus)[0], 1.4, quad(0.2, 0.2),   12., quad(12., 10.));
    m.setFreeAmplitude(m.isp()->freeAmplitudes(f_0_1370, piPlus)[0], 1.3, quad(0.4, 0.2),  -21., quad(15., 14.));
    m.setFreeAmplitude(m.isp()->freeAmplitudes(f_0_1500, piPlus)[0], 1.1, quad(0.3, 0.2),  -44., quad(13., 16.));
    m.setFreeAmplitude(m.isp()->freeAmplitudes(sigma,    piPlus)[0], 3.7, quad(0.3, 0.2),   -3., quad(4., 2.));

    // add shape parameters
    m.addParameter("f_0_980_mass", f_0_980->mass(), 0.953 - 3 * 0.02, 0.953 + 3 * 0.02);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.953, 0.02));
    m.addParameter("f_0_980_coupling", f_0_980_flatte->channels()[0].Coupling, 0.329 - 3 * 0.096, 0.329 + 3 * 0.096);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.329, 0.096));
    m.addParameter("f_0_1370_mass", f_0_1370->mass(), 1.259 - 3 * 0.055, 1.259 + 3 * 0.055);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.259, 0.055));
    m.addParameter("f_0_1370_width", std::dynamic_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->width(),
                   0.298 - 3 * 0.021, 0.298 + 3 * 0.021);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.298, 0.021));
    m.addParameter("sigma_mass", std::dynamic_pointer_cast<yap::PoleMass>(sigma->massShape())->mass(),
                   std::complex<double>(0.466, -0.223) - 3. * std::complex<double>(0.018, 0.028),
                   std::complex<double>(0.466, -0.223) + 3. * std::complex<double>(0.018, 0.028));
    m.GetParameters()[m.GetParameters().Size() - 2].SetPrior(new BCGaussianPrior(0.466, 0.018));
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(-0.223, 0.028));

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(m.axes(), m.isp(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()),
                              m.isp()->mass()->value(), m.axes(), m2r, g,
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
