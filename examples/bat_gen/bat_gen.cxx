// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_gen.h"

#include <DataSet.h>
#include <ParticleCombination.h>

#include <limits>

// -----------------------
bat_gen::bat_gen(std::string name, std::unique_ptr<yap::Model> M,
                 std::vector<std::vector<unsigned> > pcs)
    : BCModel(name),
      M_(std::move(M))
{
    if (!M_ or !M_->consistent())
        throw std::exception();

    MassAxes_ = M_->massAxes(pcs);

    for (auto& pc : MassAxes_) {
        std::string axis_label = "m2_" + indices_string(*pc).substr(1, 2);
        auto mrange = M_->massRange(pc);
        AddParameter(axis_label, pow(mrange[0], 2), pow(mrange[1], 2));
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << pow(mrange[0], 2) << ", " << pow(mrange[1], 2) << "]"
                  << std::endl;
    }
}

// ---------------------------------------------------------
double bat_gen::LogLikelihood(const std::vector<double>& parameters)
{
    unsigned c = GetCurrentChain();
    return M_->logOfSquaredAmplitude(M_->dataSet()[c], c);
}

// ---------------------------------------------------------
double bat_gen::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = M_->calculateFourMomenta(MassAxes_, parameters);

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    unsigned c = GetCurrentChain();
    M_->setFinalStateMomenta(M_->dataSet()[c], P, c);
    return 0;
}
