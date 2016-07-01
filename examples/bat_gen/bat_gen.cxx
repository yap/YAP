// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_gen.h"

#include <CalculationStatus.h>
#include <DataSet.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FreeAmplitude.h>
#include <MassRange.h>
#include <ParticleCombination.h>
#include <VariableStatus.h>

#include <complex>
#include <limits>

// -----------------------
bat_gen::bat_gen(std::string name, std::unique_ptr<yap::Model> M, std::vector<std::vector<unsigned> > pcs)
    : BCModel(name),
      Model_(std::move(M))
{
    if (!Model_ or !Model_->consistent())
        throw std::exception();

    auto isps = full_final_state_isp(*Model_);
    if (isps.empty())
        throw yap::exceptions::Exception("no full-final-state initial-state particle in model", "bat_gen::bat_gen");

    ISP_ = isps[0];

    for (auto& kv : Model_->initialStateParticles()) {
        std::cout << "Initial state particle " << to_string(*kv.first) << " with beta^2 = " << kv.second->value() << ":\n";

        auto freeAmps = freeAmplitudes(kv.first->decayTrees());

        std::cout << std::endl;
        for (const auto& fa : freeAmps)
            if (fa->variableStatus() != yap::VariableStatus::fixed)
                std::cout << to_string(*fa) << "  =  (" << abs(fa->value()) << ", " << yap::deg(arg(fa->value())) << " deg)" << std::endl;
        std::cout << std::endl;
    }

    MassAxes_ = Model_->massAxes(pcs);
    auto m2r = yap::squared(yap::mass_range(MassAxes_, ISP_, Model_->finalStateParticles()));

    for (size_t i = 0; i < MassAxes_.size(); ++i) {
        std::string axis_label = "m2_" + indices_string(*MassAxes_[i]).substr(1, 2);
        AddParameter(axis_label, m2r[i][0], m2r[i][1], axis_label, "[GeV]");
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << m2r[i][0] << ", " << m2r[i][1] << "]"
                  << std::endl;
    }

}

// ---------------------------------------------------------
void bat_gen::MCMCUserInitialize()
{
    Data_.assign(GetNChains(), Model_->createDataSet(1));
    LikelihoodCalls_.assign(GetNChains(), 0);
}

// ---------------------------------------------------------
double bat_gen::LogLikelihood(const std::vector<double>&)
{
    unsigned c = GetCurrentChain();
    double L = sum_of_log_intensity(*Model_, Data_[c]);
    // Model_->setParameterFlagsToUnchanged();
    ++LikelihoodCalls_[c];
    return L;
}

// // ---------------------------------------------------------
// void bat_gen::CalculateObservables(const std::vector<double>& )
// {
//     unsigned c = GetCurrentChain();
//     auto P = Model_->fourMomenta()->finalStateMomenta(Data_[c][0]);
//     std::vector<double> E(P.size(), 0);
//     for (size_t i = 0; i < P.size(); ++i)
//         E[i] = P[i][0] - abs(P[i]);
//     auto Esum = std::accumulate(E.begin(), E.end(), 0.);
//     for (size_t i = 0; i < E.size(); ++i)
//         GetObservable(i) = E[i] / Esum;
// }

// ---------------------------------------------------------
double bat_gen::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = Model_->calculateFourMomenta(MassAxes_, parameters, ISP_->mass()->value());

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    unsigned c = GetCurrentChain();
    Data_[c].setAll(yap::VariableStatus::changed);
    Data_[c].setAll(yap::CalculationStatus::uncalculated);
    Model_->setFinalStateMomenta(Data_[c][0], P, Data_[c]);
    return 0;
}
