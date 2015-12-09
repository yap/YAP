// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "d3pi.h"
#include "dkkpi.h"

int main()
{
    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // d3pi m("D3PI");
    dkkpi m("DKKPI");

    // set precision
    m.MCMCSetPrecision(BCEngineMCMC::kMedium);
    m.MCMCSetNChains(1);
    m.MCMCSetNIterationsPreRunCheck(1000);
    m.MCMCSetPreRunCheckClear(0);
    m.MCMCSetNIterationsPreRunMax(100000);
    m.MCMCSetNIterationsPreRunMin(5000);
    m.MCMCSetNIterationsRun(1000000);
    m.MCMCSetMultivariateProposalFunctionCovarianceUpdatesMinimum(10);

    m.GetObservables().FillHistograms(true, true);

    m.MCMCSetFlagInitialPosition(BCEngineMCMC::kMCMCInitRandomUniform);

    m.MCMCSetMultivariateProposalFunction(true);
    m.MCMCSetMultivariateProposalFunctionCovarianceUpdateLambda(0.5);
    m.MCMCSetMinimumEfficiency(0.15);
    m.MCMCSetMaximumEfficiency(0.35);
    m.MCMCSetRValueParametersCriterion(1.25);

    m.MCMCSetCorrectRValueForSamplingVariability(true);

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
