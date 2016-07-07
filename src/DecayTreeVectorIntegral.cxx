#include "DecayTreeVectorIntegral.h"

#include "Constants.h"
#include "DecayTree.h"

#include <algorithm>
#include <numeric>

namespace yap {

//-------------------------
DecayTreeVectorIntegral::DecayTreeVectorIntegral(const DecayTreeVector& dtv)
    : DecayTrees_(dtv)
{
    // initialize diagonal elements
    for (const auto& dt : DecayTrees_)
        Diagonals_.emplace(dt, DiagonalIntegralMap::mapped_type());
    // initialize off-diagonal elements
    for (size_t i = 0; i < DecayTrees_.size(); ++i)
        for (size_t j = i + 1; j < DecayTrees_.size(); ++j)
            OffDiagonals_.emplace(OffDiagonalIntegralMap::key_type({DecayTrees_[i], DecayTrees_[j]}),
                                  OffDiagonalIntegralMap::mapped_type());
}

//-------------------------
const Model* DecayTreeVectorIntegral::model() const
{
    return (!DecayTrees_.empty() and DecayTrees_[0]) ? DecayTrees_[0]->model() : nullptr;
}

//-------------------------
const std::vector<double> fit_fractions(const DecayTreeVectorIntegral& MI)
{
    double I = MI.integral().value;
    std::vector<double> ff;
    ff.reserve(MI.decayTrees().size());
    std::transform(MI.decayTrees().begin(), MI.decayTrees().end(), std::back_inserter(ff),
    [&](const DecayTreeVector::value_type & dt) {return integral(*MI.diagonals().find(dt)).value / I;});
    return ff;
}

//-------------------------
const std::vector<std::vector<std::complex<double> > > cached_integrals(const DecayTreeVectorIntegral& MI)
{
    std::vector<std::vector<std::complex<double> > > I(MI.decayTrees().size(), std::vector<std::complex<double> >(MI.decayTrees().size(), Complex_0));
    for (size_t i = 0; i < MI.decayTrees().size(); ++i) {
        I[i][i] = MI.diagonals().at(MI.decayTrees()[i]).value;
        for (size_t j = i + 1; j < MI.decayTrees().size(); ++j)
            I[j][i] = conj(I[i][j] = MI.offDiagonals().at({MI.decayTrees()[i], MI.decayTrees()[j]}).value);
    }
    return I;
}

//-------------------------
const std::vector<std::vector<std::complex<double> > > integrals(const DecayTreeVectorIntegral& MI)
{
    std::vector<std::vector<std::complex<double> > > I(MI.decayTrees().size(), std::vector<std::complex<double> >(MI.decayTrees().size(), Complex_0));
    for (size_t i = 0; i < MI.decayTrees().size(); ++i) {
        I[i][i] = integral(*MI.diagonals().find(MI.decayTrees()[i])).value;
        for (size_t j = i + 1; j < MI.decayTrees().size(); ++j)
            I[j][i] = conj(I[i][j] = integral(*MI.offDiagonals().find({MI.decayTrees()[i], MI.decayTrees()[j]})).value);
    }
    return I;
}

//-------------------------
const RealIntegralElement integral(const DiagonalIntegralMap::value_type& a_A2)
{
    return RealIntegralElement(norm(a_A2.first->dataIndependentAmplitude()) * a_A2.second.value);
}

//-------------------------
const RealIntegralElement integral(const OffDiagonalIntegralMap::value_type& aa_AA)
{
    return RealIntegralElement(real(2. * conj(aa_AA.first[0]->dataIndependentAmplitude())
                                    * aa_AA.first[1]->dataIndependentAmplitude()
                                    * aa_AA.second.value));
}

//-------------------------
const RealIntegralElement operator+(const RealIntegralElement& d, const DiagonalIntegralMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const RealIntegralElement operator+(const RealIntegralElement& d, const OffDiagonalIntegralMap::value_type& v)
{ return d + integral(v); }

//-------------------------
const IntegralElement<double> DecayTreeVectorIntegral::integral() const
{
    return std::accumulate(Diagonals_.begin(), Diagonals_.end(), RealIntegralElement())
           + std::accumulate(OffDiagonals_.begin(), OffDiagonals_.end(), RealIntegralElement());
}

}
