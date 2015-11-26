// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKKPI__H
#define __BAT__DKKPI__H

#include <BAT/BCModel.h>

#include <InitialStateParticle.h>
#include <ParticleCombination.h>

#include <memory>
#include <string>
#include <vector>

// This is a dkkpi header file.
// Model source code is located in file yap_test/dkkpi.cxx

// ---------------------------------------------------------
class dkkpi : public BCModel
{

public:

    // Constructor
    dkkpi(std::string name);

    // Destructor
    ~dkkpi();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double> & parameters);

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double> & parameters);

protected:
    yap::ParticleCombinationVector DalitzAxes_;
    std::shared_ptr<yap::InitialStateParticle> D_;

    double m2_P;
    double m2_a;
    double m2_b;
    double m2_c;

};
// ---------------------------------------------------------

#endif
