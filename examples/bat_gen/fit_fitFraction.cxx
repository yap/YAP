#include "fit_fitFraction.h"

#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <DecayTreeVectorIntegral.h>
#include <FreeAmplitude.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <Spin.h>
#include <VariableStatus.h>

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>

//-------------------------
fit_fitFraction::fit_fitFraction(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_fit(name, std::move(M), pcs)
{
    AddParameter("N_0", 1.e-3, 50.);
    AddParameter("N_1", 1.e-3, 1.e3);
    AddParameter("N_2", 1.e-3, 2);

    SetPriorConstantAll();

    FitFractions_.assign(CalculatedFitFractions_[0].size(), { -1, 0});
}

//-------------------------
void fit_fitFraction::setFitFraction(std::shared_ptr<yap::DecayTree> dt, double mean, double sigma)
{
    auto it = std::find(DecayTrees_.begin(), DecayTrees_.end(), dt);
    if (it == DecayTrees_.end())
        throw yap::exceptions::Exception("could not find decay tree", "fit_fitFraction::setFirFraction");
    FitFractions_[it - DecayTrees_.begin()] = {mean, sigma};
}

//-------------------------
double fit_fitFraction::LogLikelihood(const std::vector<double>& p)
{
    // apply normalization change
    auto P = p;
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        P[i * 2] *= p[FreeAmplitudes_.size() * 2 + FreeAmplitudes_[i]->spinAmplitude()->L()];

    setParameters(P);
    unsigned c = GetCurrentChain();

    double L = 0;
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        if (FitFractions_[i][0] > 0)
            L += BCMath::LogGaus(CalculatedFitFractions_[c][i].value(), FitFractions_[i][0], FitFractions_[i][1]);

    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();

    return L;
}

