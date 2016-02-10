// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BAT_GEN__H
#define __BAT__BAT_GEN__H

#include <BAT/BCModel.h>

#include <MassAxes.h>
#include <Model.h>
#include <ParticleFactory.h>
#include <SpinAmplitudeCache.h>

#include <memory>
#include <string>
#include <vector>

// This is a bat_gen header file.
// Model source code is located in file yap_test/bat_gen.cxx

// ---------------------------------------------------------
class bat_gen : public BCModel
{

public:

    // Constructor
    bat_gen(std::string name, std::unique_ptr<yap::Model> M,
            std::vector<std::vector<unsigned> > pcs);

    void initialize(unsigned n)
    { M_->initializeForMonteCarloGeneration(n); }

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& parameters);

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double>& parameters);

protected:
    yap::MassAxes MassAxes_;
    std::shared_ptr<yap::Model> M_;

};
// ---------------------------------------------------------

#endif
