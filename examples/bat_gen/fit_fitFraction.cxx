#include "fit_fitFraction.h"

#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <DecayTreeVectorIntegral.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <Spin.h>
#include <VariableStatus.h>

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>

//-------------------------
fit_fitFraction::fit_fitFraction(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_fit(name, std::move(M), pcs),
      FitFractions_(decayTrees().size(), { -1, 0})
{
    AddParameter("N_0", 1.e-3, 50.);
    AddParameter("N_1", 1.e-3, 1.e3);
    AddParameter("N_2", 1.e-3, 2);

    for (const auto& fa : free_amplitudes(decayTrees())) {
        if (fa->variableStatus() == yap::VariableStatus::fixed)
            continue;
        AddParameter("amp(" + to_string(*fa->decayChannel()) + " M = " + yap::spin_to_string(fa->twoM()) + ")",
                     0, 2 * norm(fa->value()));
        AddParameter("phase(" + to_string(*fa->decayChannel()) + " M = " + yap::spin_to_string(fa->twoM()) + ")",
                     -180, 180);
        FreeAmplitudes_.push_back(fa);
    }

    SetPriorConstantAll();

    for (const auto& dt : decayTrees())
        AddObservable("fit_frac(" + to_string(*dt->freeAmplitude()->decayChannel()) + " M = " + yap::spin_to_string(dt->freeAmplitude()->twoM()) + ")", 0, 1.1);
}

//-------------------------
const yap::DecayTreeVector& fit_fitFraction::decayTrees() const
{ return isp()->decayTrees().begin()->second; }

//-------------------------
const yap::DecayTreeVectorIntegral& fit_fitFraction::decayTreesIntegral() const
{ return Integral_.integrals().begin()->second; }

//-------------------------
void fit_fitFraction::setFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A, double amp, double amp_s, double phase, double phase_s)
{
    auto it = std::find(FreeAmplitudes_.begin(), FreeAmplitudes_.end(), A);

    if (it == FreeAmplitudes_.end())
        throw yap::exceptions::Exception("FreeAmplitude not found", "addParameter");

    size_t i = 3 + (it - FreeAmplitudes_.begin()) * 2;

    // add amp
    if (amp_s > 0) {
        GetParameters()[i].SetLimits(amp - 3. * amp_s, amp + 3. * amp_s);
        GetParameters()[i].SetPrior(new BCGaussianPrior(amp, amp_s));
    } else {
        GetParameters()[i].SetLimits(0, 2 * amp);
        GetParameters()[i].SetPriorConstant();
        GetParameters()[i].Fix(amp);
    }
    if (phase_s > 0) {
        GetParameters()[i + 1].SetLimits(phase - 3. * phase_s, phase + 3. * phase_s);
        GetParameters()[i + 1].SetPrior(new BCGaussianPrior(phase, phase_s));
    } else {
        GetParameters()[i + 1].SetLimits(phase - 10, phase + 10);
        GetParameters()[i + 1].SetPriorConstant();
        GetParameters()[i + 1].Fix(phase);
    }
}

//-------------------------
void fit_fitFraction::MCMCUserInitialize()
{
    bat_fit::MCMCUserInitialize();
    CalculatedFitFractions_.assign(GetNChains(), yap::RealIntegralElementVector());
}

//-------------------------
void fit_fitFraction::setFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double stat, double syst)
{
    auto it = std::find(decayTrees().begin(), decayTrees().end(), dt);
    if (it == decayTrees().end())
        throw yap::exceptions::Exception("could not find decay tree", "fit_fitFraction::setFirFraction");
    FitFractions_[it - decayTrees().begin()] = {mean, sqrt(stat * stat + syst * syst)};
}

//-------------------------
double fit_fitFraction::LogLikelihood(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = p[FreeAmplitudes_[i]->spinAmplitude()->L()] * std::polar(p[3 + i * 2], yap::rad(p[3 + i * 2 + 1]));

    yap::set_values(Parameters_.begin(), Parameters_.end(), p.begin() + 3 + FreeAmplitudes_.size() * 2, p.end());

    Integrator_(Integral_, IntegralPartitions_);

    unsigned c = GetCurrentChain();
    CalculatedFitFractions_[c] = fit_fractions(decayTreesIntegral());

    double L = 0;
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        if (FitFractions_[i][0] > 0)
            L += BCMath::LogGaus(CalculatedFitFractions_[c][i].value(), FitFractions_[i][0], FitFractions_[i][1]);

    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();

    return L;
}

//-------------------------
void fit_fitFraction::CalculateObservables(const std::vector<double>& p)
{
    unsigned c = GetCurrentChain();
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        GetObservables()[i] = CalculatedFitFractions_[c][i].value();
}
