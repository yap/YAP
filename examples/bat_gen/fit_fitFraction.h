#ifndef __BAT__fit_frac__H
#define __BAT__fit_frac__H

#include <fwd/DecayTree.h>

#include "bat_fit.h"

#include <array>
#include <memory>
#include <vector>

class fit_fitFraction : public bat_fit
{

public:
    /// constructor
    fit_fitFraction(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {});

    /// log likelihood
    double LogLikelihood(const std::vector<double>&) override;

    /// set fit fraction
    void setFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double sigma);

private:

    /// Fit fraction data to fit to
    std::vector<std::array<double, 2> > FitFractions_;

};

#endif
