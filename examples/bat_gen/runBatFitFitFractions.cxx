// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d3pi.h"
#include "models/d4pi.h"
#include "tools.h"

#include <DecayTree.h>
#include <Filters.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <Model.h>
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

int main()
{
    yap::plainLogs(el::Level::Info);

    // create bat_fit object
    auto m = d3pi_fit_fitFractions("D3PI_frac_fit", yap_model<yap::ZemachFormalism>());
    double D_mass = 1.86961;
    //auto m = d4pi_fit_fitFraction();
    //double D_mass = 1.8648400;

    FLOG(INFO) << "Created model";

    m.GetParameter("N_1").Fix(1);

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), D_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
    m.setNIntegrationPoints(1e6, 1e5);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e7);
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

    // for (size_t i = 0; i < D->decayTrees().at(0).size(); ++i)
    //     std::cout << i << " = " << to_string(*D->decayTrees().at(0)[i]) << std::endl;

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    // m.PrintKnowledgeUpdatePlots("output/" + m.GetSafeName() + "_update.pdf", 2, 2, true);

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
