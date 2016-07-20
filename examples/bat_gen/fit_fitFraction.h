#ifndef __BAT__fit_frac__H
#define __BAT__fit_frac__H

#include <DecayTree.h>
#include <FreeAmplitude.h>

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

    void CalculateObservables(const std::vector<double>&) override;

    using fit_fraction_map = std::map<std::shared_ptr<yap::DecayTree>, std::array<double, 2> >;

    /// add a FreeAmplitude
    void addFreeAmplitude(std::string name, std::shared_ptr<yap::FreeAmplitude> A, double amp, double amp_s, double phase, double phase_s);
    void addFreeAmplitude(std::string name, std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
    { addFreeAmplitude(name, A, amp, -1, phase, -1); }

    /// add fit fraction
    void addFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double stat, double syst);
    
private:

    yap::DecayTreeVector DecayTrees_;

    fit_fraction_map FitFractions_;
  
    yap::FreeAmplitudeVector FreeAmplitudes_;

};

#endif
