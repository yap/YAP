// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include <BAT/BCModel.h>

#include <InitialStateParticle.h>
#include <MassAxes.h>

#include <memory>
#include <string>
#include <vector>

// This is a d3pi header file.
// Model source code is located in file yap_test/d3pi.cxx

// ---------------------------------------------------------
class d3pi : public BCModel
{

public:

    // Constructor
    d3pi(std::string name);

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& parameters);

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double>& parameters);

protected:
    yap::MassAxes MassAxes_;
    std::shared_ptr<yap::InitialStateParticle> D_;

};
// ---------------------------------------------------------

#endif
