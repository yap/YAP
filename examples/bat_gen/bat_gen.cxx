// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_gen.h"

#include <DataSet.h>
#include <ParticleCombination.h>
#include <CalculationStatus.h>

#include <limits>

// -----------------------
bat_gen::bat_gen(std::string name, std::unique_ptr<yap::Model> M, std::vector<std::vector<unsigned> > pcs)
    : BCModel(name),
      Model_(std::move(M))
{
    if (!Model_ or !Model_->consistent())
        throw std::exception();

    MassAxes_ = Model_->massAxes(pcs);

    for (auto& pc : MassAxes_) {
        std::string axis_label = "m2_" + indices_string(*pc).substr(1, 2);
        auto mrange = Model_->massRange(pc);
        AddParameter(axis_label, pow(mrange[0], 2), pow(mrange[1], 2));
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << pow(mrange[0], 2) << ", " << pow(mrange[1], 2) << "]"
                  << std::endl;
    }
}

// ---------------------------------------------------------
void bat_gen::MCMCUserInitialize()
{
    Data_.clear();
    Data_.assign(GetNChains(), Model_->dataSet());
    for (auto& d : Data_)
        d.addEmptyPoint();

    // // create a data partition for each chain
    // Partitions_ = yap::DataPartitionBlock::createBySize(Data_, 1);
}

// ---------------------------------------------------------
double bat_gen::LogLikelihood(const std::vector<double>& parameters)
{
    unsigned c = GetCurrentChain();
    return Model_->sumOfLogsOfSquaredAmplitudes(Data_[c]);
    // return Model_->partialSumOfLogsOfSquaredAmplitudes(Partitions_[c].get(), Data_);
}

// ---------------------------------------------------------
double bat_gen::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = Model_->calculateFourMomenta(MassAxes_, parameters);

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    unsigned c = GetCurrentChain();
    Data_[c][0].setFinalStateMomenta(P);
    return 0;
}
