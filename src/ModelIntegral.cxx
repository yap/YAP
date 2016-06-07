#include "ModelIntegral.h"

#include "DecayTree.h"

#include <numeric>

namespace yap {

//-------------------------
ModelIntegral::ModelIntegral(const DecayTreeVector& dtv)
    : DecayTrees_(dtv)
{
    // initialize diagonal elements
    for (const auto& dt : DecayTrees_)
        Diagonals_.emplace(dt, 0);
    // initialize off-diagonal elements
    for (size_t i = 0; i < DecayTrees_.size(); ++i)
        for (size_t j = i + 1; j < DecayTrees_.size(); ++j)
            OffDiagonals_.emplace(OffDiagonalMap::key_type({DecayTrees_[i], DecayTrees_[j]}), 0.);
}

//-------------------------
const double integral(const ModelIntegral::DiagonalMap::value_type& a_A2)
{
    return norm(a_A2.first->dataIndependentAmplitude()) * a_A2.second;
}

//-------------------------
const double integral(const ModelIntegral::OffDiagonalMap::value_type& aa_AA)
{
    return 2. * real(aa_AA.first[0]->dataIndependentAmplitude() *
                     aa_AA.first[1]->dataIndependentAmplitude() *
                     aa_AA.second);
}

//-------------------------
const double operator+(const double& d, const ModelIntegral::DiagonalMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const double operator+(const double& d, const ModelIntegral::OffDiagonalMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const double ModelIntegral::integral() const
{
    return std::accumulate(Diagonals_.begin(), Diagonals_.end(), 0.)
           + std::accumulate(OffDiagonals_.begin(), OffDiagonals_.end(), 0.);
}

}
