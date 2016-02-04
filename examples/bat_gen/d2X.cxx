// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "d2X.h"

#include <BAT/BCMath.h>

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"

#include <complex>
#include <limits>

// ---------------------------------------------------------
d2X::d2X(std::string name)
    : BCModel(name),
      Factory_((std::string)::getenv("YAPDIR") + "/evt.pdl")
{
}

//-------------------------
void d2X::initialize(unsigned n)
{
    if (!D_)
        throw std::exception();

    D_->initializeForMonteCarloGeneration(n);
    std::cout << "number of data partitions = " << D_->dataPartitions().size() << std::endl;

    MassAxes_ = D_->getMassAxes({{0, 1}, {1, 2}});

    for (auto& pc : MassAxes_) {
        std::string axis_label = "m2_" + indices_string(*pc).substr(1, 2);
        auto mrange = D_->getMassRange(pc);
        AddParameter(axis_label, pow(mrange[0], 2), pow(mrange[1], 2));
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << pow(mrange[0], 2) << ", " << pow(mrange[1], 2) << "]"
                  << std::endl;
    }
}

// ---------------------------------------------------------
double d2X::LogLikelihood(const std::vector<double>& parameters)
{
    unsigned c = GetCurrentChain();
    return D_->logOfSquaredAmplitude(D_->dataSet()[c], c);
}

// ---------------------------------------------------------
double d2X::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = D_->calculateFourMomenta(MassAxes_, parameters);

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    unsigned c = GetCurrentChain();
    D_->setFinalStateFourMomenta(D_->dataSet()[c], P, c);
    return 0;
}
