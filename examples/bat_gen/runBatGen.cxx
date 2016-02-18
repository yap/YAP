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
#include <make_unique.h>
#include <ZemachFormalism.h>

#include "bat_gen.h"
#include "d3pi.h"
#include "dkkpi.h"

int main()
{
    yap::plainLogs(el::Level::Info);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // bat_gen m("D3PI", std::move(d3pi(std::make_unique<yap::ZemachFormalism>())), {{0, 1}, {1, 2}});
    bat_gen m("DKKPI", std::move(dkkpi(std::make_unique<yap::ZemachFormalism>())), {{0, 1}, {1, 2}});
    // bat_gen m("DKKPI", std::move(dkkpi(std::make_unique<yap::HelicityFormalism>())), {{0, 1}, {1, 2}});

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNChains(4);

    BCLog::OutSummary("Initializing for MC Generation");
    m.initialize(m.GetNChains());

    m.SetNIterationsRun(static_cast<int>(5e5 / m.GetNChains()));

    m.GetObservables().FillHistograms(true, true);

    BCLog::OutSummary("Test model created");

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full par. space
    // m.Normalize();

    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    // m.FindMode(m.GetGlobalMode());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary plots
    // m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    // m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
