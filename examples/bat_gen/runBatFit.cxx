// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
//#include "hist.h"
#include "models/d3pi.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"

#include <HelicityFormalism.h>
#include <MassRange.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <make_unique.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

template <typename T>
std::unique_ptr<yap::Model> fit_model(bool include_phsp_factors = true)
{ return std::make_unique<yap::Model>(std::make_unique<T>(), include_phsp_factors); }

int main()
{
    yap::plainLogs(el::Level::Info);

    // open file
    // std::string model_name = "D3PI";
    //std::string model_name = "DKKPI";
    std::string model_name = "D4PI";

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

    // create model
    // auto m = d3pi_fit(model_name + "_fit", fit_model<yap::ZemachFormalism>(false), find_mass_axes(*t_pars));
    // auto m = dkkpi_fit(model_name + "_fit", fit_model<yap::HelicityFormalism>(false), find_mass_axes(*t_pars));
    auto m = d4pi_fit(model_name + "_fit", find_mass_axes(*t_pars));

    //double D_mass = 1.86961;
    double D_mass = 1.8648400;

    // load fit data and partition it
    load_data(m.fitData(), *m.model(), m.axes(), D_mass, *t_mcmc, 20000, 45);
    // m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), 2);

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 40000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), D_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max()));
    LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";

    // partition integration data
    // m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 2);

    // TH2D* h2_fit_data = hist2(*m.model()->fourMomenta(), m.axes(), m2r, m.fitData());
    // TH2D* h2_int_data = hist2(*m.model()->fourMomenta(), m.axes(), m2r, m.integralData());

    // TH1D* h1_fit_data = hist1(*m.model()->fourMomenta(), m.axes()[0], m2r[0], m.fitData());
    // TH1D* h1_int_data = hist1(*m.model()->fourMomenta(), m.axes()[0], m2r[0], m.integralData());

    // TCanvas* C = new TCanvas("C","C");
    // C->Divide(1,2);
    // C->cd(1);
    // h2_fit_data->Draw("colz");
    // C->cd(2);
    // h2_int_data->Draw("colz");
    // C->Print("data.pdf");

    // return 0;

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

    m.FindMode(m.GetBestFitParameters());

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    // m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    // m.PrintKnowledgeUpdatePlots("output/" + m.GetSafeName() + "_update.pdf", 2, 2, false);//true);

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
