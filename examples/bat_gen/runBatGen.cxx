// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <HelicityFormalism.h>
#include <make_unique.h>
#include <ZemachFormalism.h>

#include "bat_gen.h"
#include "d3pi.h"
#include "dkkpi.h"
#include "D_K0pi0pi0.h"

#include <chrono>
#include <ratio>

int main()
{
    yap::plainLogs(el::Level::Info);

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // bat_gen m("D3PI", std::move(d3pi(std::make_unique<yap::ZemachFormalism>())), {{0, 1}, {1, 2}});
    bat_gen m("DKSPIPI", std::move(D_K0pi0pi0(std::make_unique<yap::ZemachFormalism>())), {{0, 1}, {1, 2}});
    // bat_gen m("DKKPI", std::move(dkkpi(std::make_unique<yap::ZemachFormalism>())), {{0, 1}, {1, 2}});
    // bat_gen m("DKKPI", std::move(dkkpi(std::make_unique<yap::HelicityFormalism>())), {{0, 1}, {1, 2}});

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNChains(4);
    m.SetMinimumEfficiency(0.85);
    m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(1e6 / m.GetNChains()));

    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    // m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
