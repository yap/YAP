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
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNChains(1);
    m.SetNIterationsPreRunCheck(1000);
    m.SetPreRunCheckClear(0);
    m.SetNIterationsPreRunMax(100000);
    m.SetNIterationsPreRunMin(5000);
    m.SetNIterationsRun(1000000);
    // \todo has been renamed in latest BAT version
    // m.SetMultivariateProposalFunctionCovarianceUpdatesMinimum(10);

    m.GetObservables().FillHistograms(true, true);

    // \todo has been renamed in latest BAT version
    //m.SetFlagInitialPosition(BCEngineMCMC::kMCMCInitRandomUniform);

    // \todo has been renamed in latest BAT version
    //m.SetMultivariateProposalFunction(true);
    // \todo has been renamed in latest BAT version
    //m.SetMultivariateProposalFunctionCovarianceUpdateLambda(0.5);
    m.SetMinimumEfficiency(0.15);
    m.SetMaximumEfficiency(0.35);
    m.SetRValueParametersCriterion(1.25);

    m.SetCorrectRValueForSamplingVariability(true);

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
