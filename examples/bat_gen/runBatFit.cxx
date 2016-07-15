// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
#include "models/d3pi.h"

#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <PHSP.h>
#include <ZemachFormalism.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    // open file
    std::string model_name = "D3PI";
    auto file = TFile::Open(("output/" + model_name + "_mcmc.root").data(), "READ");
    if (file->IsZombie())
        throw yap::exceptions::Exception("could not open file", "main");

    TTree* t_mcmc = nullptr;
    file->GetObject((model_name + "_mcmc").data(), t_mcmc);
    if (!t_mcmc)
        throw yap::exceptions::Exception("could not retrieve mcmc tree", "main");

    TTree* t_pars = nullptr;
    file->GetObject((model_name + "_parameters").data(), t_pars);
    if (!t_pars)
        throw yap::exceptions::Exception("could not retrieve mcmc tree", "main");
    
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

    M->addInitialStateParticle(D);

    // create bat_fit object
    bat_fit m("D3PI_fit", std::move(M), *t_pars);

    // set parameters of fit
    auto amp_range = std::complex<double>(5., 5.);
    D->freeAmplitudes(rho, piPlus)[0]->variableStatus() = yap::VariableStatus::fixed;
    m.addParameter("f_0_980",  m.isp()->freeAmplitudes(f_0_980,  piPlus)[0], -amp_range, amp_range);
    m.addParameter("f_2",      m.isp()->freeAmplitudes(f_2,      piPlus)[0], -amp_range, amp_range);
    m.addParameter("f_0_1370", m.isp()->freeAmplitudes(f_0_1370, piPlus)[0], -amp_range, amp_range);
    m.addParameter("f_0_1500", m.isp()->freeAmplitudes(f_0_1500, piPlus)[0], -amp_range, amp_range);
    m.addParameter("sigma",    m.isp()->freeAmplitudes(sigma,    piPlus)[0], -amp_range, amp_range);

    m.GetParameters().SetPriorConstantAll();

    // load fit data
    m.loadData(m.fitData(), *t_mcmc, 10000, 1);

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(m.axes(), m.isp(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()),
                              m.isp()->mass()->value(), m.axes(), m2r, g,
                              std::numeric_limits<unsigned>::max()));
    LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(1e5 / m.GetNChains()));

    m.WriteMarkovChain("output/" + m.GetSafeName() + "_mcmc.root", "RECREATE");

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

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
