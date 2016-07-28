#ifndef __BAT__fit_frac__H
#define __BAT__fit_frac__H

#include <fwd/DecayTree.h>
#include <fwd/DecayTreeVectorIntegral.h>
#include <fwd/IntegralElement.h>

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
    void setFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A, double amp, double amp_s, double phase, double phase_s);
    void setFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
    { setFreeAmplitude(A, amp, -1, phase, -1); }

    /// set fit fraction
    void setFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double stat, double syst);

    yap::RealIntegralElementVector setParameters(const std::vector<double>& p);

    const yap::DecayTreeVector& decayTrees() const;

    const yap::DecayTreeVectorIntegral& decayTreesIntegral() const;

private:

    std::vector<std::array<double, 2> > FitFractions_;

    yap::FreeAmplitudeVector FreeAmplitudes_;

};

#endif
