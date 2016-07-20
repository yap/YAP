#include "fit_fitFraction.h"

#include <DecayingParticle.h>
#include <DecayTreeVectorIntegral.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>

//-------------------------
fit_fitFraction::fit_fitFraction(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_fit(name, std::move(M), pcs),
      DecayTrees_(isp()->decayTrees().begin()->second)
{
    AddParameter("N_0", 1.e-3, 1.e3);
    AddParameter("N_1", 1.e-3, 1.e3);
    AddParameter("N_2", 1.e-3, 1.e3);
    SetPriorConstantAll();
}

//-------------------------
void fit_fitFraction::addFreeAmplitude(std::string name, std::shared_ptr<yap::FreeAmplitude> A, double amp, double amp_s, double phase, double phase_s)
{
    if (std::find(FreeAmplitudes_.begin(), FreeAmplitudes_.end(), A) != FreeAmplitudes_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "addParameter");
    FreeAmplitudes_.push_back(A);
    // add amp
    if (amp_s > 0) {
        AddParameter(name + "_amp", amp - 3 * amp_s, amp + 3 * amp_s);
        GetParameters().Back().SetPrior(new BCGaussianPrior(amp, amp_s));
        // GetParameters().Back().Fix(amp);
    } else {
        AddParameter(name + "_amp", amp - 1, amp + 1);
        GetParameters().Back().SetPriorConstant();
        GetParameters().Back().Fix(amp);
    }
    if (phase_s > 0) {
        AddParameter(name + "_phs", phase - 3 * phase_s, phase + 3 * phase_s);
        GetParameters().Back().SetPrior(new BCGaussianPrior(phase, phase_s));
        // GetParameters().Back().Fix(phase);
    } else {
        AddParameter(name + "_phs", phase - 0.5, phase + 0.5);
        GetParameters().Back().SetPrior(new BCGaussianPrior(phase, phase_s));
        GetParameters().Back().Fix(phase);
    }
}

//-------------------------
void fit_fitFraction::addFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double stat, double syst)
{
    FitFractions_[dt] = {mean, sqrt(stat * stat + syst * syst)};
    AddObservable("ff_" + std::to_string(GetObservables().Size()), 0, 2);
}

//-------------------------
double fit_fitFraction::LogLikelihood(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = p[FreeAmplitudes_[i]->spinAmplitude()->L()] * std::polar(p[3 + i * 2 + 0], p[3 + i * 2 + 1]);
    
    Integrator_(Integral_, IntegralPartitions_);

    auto ff = fit_fractions(Integral_.integral(DecayTrees_));
    double L = 0;
    for (size_t i = 0; i < ff.size(); ++i) {
        const auto& m_s = FitFractions_.at(DecayTrees_[i]);
        L += BCMath::LogGaus(ff[i], m_s[0], m_s[1]);
    }

    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();

    return L;
}

//-------------------------
void fit_fitFraction::CalculateObservables(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = p[FreeAmplitudes_[i]->spinAmplitude()->L()] * std::polar(p[3 + i * 2 + 0], p[3 + i * 2 + 1]);
    
    Integrator_(Integral_, IntegralPartitions_);
    
    auto ff = fit_fractions(Integral_.integral(DecayTrees_));
    for (size_t i = 0; i < ff.size(); ++i)
        GetObservables()[i] = ff[i];
}
