#include "ModelIntegral.h"

#include "DecayTree.h"

#include <algorithm>
#include <numeric>

namespace yap {

//-------------------------
ModelIntegral::ModelIntegral(const DecayTreeVector& dtv)
    : DecayTrees_(dtv)
{
    // initialize diagonal elements
    for (const auto& dt : DecayTrees_)
        Diagonals_.emplace(dt, DiagonalMap::mapped_type());
    // initialize off-diagonal elements
    for (size_t i = 0; i < DecayTrees_.size(); ++i)
        for (size_t j = i + 1; j < DecayTrees_.size(); ++j)
            OffDiagonals_.emplace(OffDiagonalMap::key_type({DecayTrees_[i], DecayTrees_[j]}),
                                  OffDiagonalMap::mapped_type());
}

//-------------------------
const Model* ModelIntegral::model() const
{
    return (!DecayTrees_.empty() and DecayTrees_[0]) ? DecayTrees_[0]->model() : nullptr;
}

//-------------------------
const std::vector<double> fitFractions(const ModelIntegral& MI)
{
    double I = MI.integral().value;
    std::vector<double> ff;
    ff.reserve(MI.decayTrees().size());
    std::transform(MI.decayTrees().begin(), MI.decayTrees().end(), std::back_inserter(ff),
    [&](const DecayTreeVector::value_type & dt) {return integral(*MI.diagonals().find(dt)).value / I;});
    return ff;
}

//-------------------------
const RealIntegralElement integral(const ModelIntegral::DiagonalMap::value_type& a_A2)
{
    return RealIntegralElement(norm(a_A2.first->dataIndependentAmplitude()) * a_A2.second.value);
}

//-------------------------
const RealIntegralElement integral(const ModelIntegral::OffDiagonalMap::value_type& aa_AA)
{
    return RealIntegralElement(real(2. * conj(aa_AA.first[0]->dataIndependentAmplitude())
                                    * aa_AA.first[1]->dataIndependentAmplitude()
                                    * aa_AA.second.value));
}

//-------------------------
const RealIntegralElement operator+(const RealIntegralElement& d, const ModelIntegral::DiagonalMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const RealIntegralElement operator+(const RealIntegralElement& d, const ModelIntegral::OffDiagonalMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const IntegralElement<double> ModelIntegral::integral() const
{
    return std::accumulate(Diagonals_.begin(), Diagonals_.end(), RealIntegralElement())
           + std::accumulate(OffDiagonals_.begin(), OffDiagonals_.end(), RealIntegralElement());
}

}
