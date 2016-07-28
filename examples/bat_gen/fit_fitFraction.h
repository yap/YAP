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

    /// calculate observables
    void CalculateObservables(const std::vector<double>&) override;

    /// set the range and prior for a FreeAmplitude
    void setFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A, double amp, double amp_s, double phase, double phase_s);

    /// fix a FreeAmplitude
    void fixFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
    { setFreeAmplitude(A, amp, -1, phase, -1); }

    /// set fit fraction
    void setFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double stat, double syst);

    /// \return decay trees of fit
    const yap::DecayTreeVector& decayTrees() const;

    /// \return decay tree integral of fit
    const yap::DecayTreeVectorIntegral& decayTreesIntegral() const;

    /// init size of CalculatedFitFractions_
    void MCMCUserInitialize();

private:

    /// Fit fraction data to fit to
    std::vector<std::array<double, 2> > FitFractions_;

    /// Calculated fit fractions (for observables)
    std::vector<yap::RealIntegralElementVector> CalculatedFitFractions_;

    /// Free amplitudes of model to set
    yap::FreeAmplitudeVector FreeAmplitudes_;

};

#endif
