// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D2X__H
#define __BAT__D2X__H

#include <BAT/BCModel.h>

#include <InitialStateParticle.h>
#include <MassAxes.h>
#include <ParticleFactory.h>

#include <memory>
#include <string>
#include <vector>

// This is a d2X header file.
// Model source code is located in file yap_test/d2X.cxx

// ---------------------------------------------------------
class d2X : public BCModel
{

public:

    // Constructor
    d2X(std::string name);

    void initialize(unsigned n);

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& parameters);

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double>& parameters);

protected:
    yap::MassAxes MassAxes_;
    std::shared_ptr<yap::InitialStateParticle> D_;
    yap::ParticleFactory Factory_;

};
// ---------------------------------------------------------

#endif
