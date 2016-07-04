// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "bat_fit.h"
#include "models/d3pi.h"

#include <logging.h>
#include <make_unique.h>
#include <ZemachFormalism.h>

#include <TFile.h>
#include <TTree.h>

#include <chrono>

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

    bat_fit m("D3PI", d3pi(std::make_unique<yap::ZemachFormalism>()),
              *t_mcmc, *t_pars,
              100, 1);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // // set precision
    // m.SetPrecision(BCEngineMCMC::kMedium);
    // m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    // m.SetNIterationsRun(static_cast<int>(1e6 / m.GetNChains()));

    // m.WriteMarkovChain("output/" + m.GetSafeName() + "_mcmc.root", "RECREATE");

    // // start timing:
    // auto start = std::chrono::steady_clock::now();

    // // run MCMC, marginalizing posterior
    // m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // // end timing
    // auto end = std::chrono::steady_clock::now();

    // m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // // timing:
    // auto diff = end - start;
    // auto ms = std::chrono::duration<double, std::micro>(diff).count();
    // auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    // BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    // BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    // BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
