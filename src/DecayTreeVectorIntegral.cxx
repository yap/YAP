#include "DecayTreeVectorIntegral.h"

#include "DecayTree.h"
#include "Exceptions.h"
#include "IntegralElement.h"

#include <algorithm>
#include <complex>
#include <numeric>

namespace yap {

//-------------------------
DecayTreeVectorIntegral::DecayTreeVectorIntegral(const DecayTreeVector& dtv)
    : DecayTrees_(dtv),
      Diagonals_(DecayTrees_.size()),
      OffDiagonals_(DecayTrees_.size() - 1)
{
    for (size_t i = 0; i < Diagonals_.size(); ++i)
        Diagonals_[i].value() = 1.;

    for (size_t i = 0; i < OffDiagonals_.size(); ++i)
        OffDiagonals_[i] = ComplexIntegralElementMatrix::value_type(DecayTrees_.size() - i - 1);
}

//-------------------------
const Model* DecayTreeVectorIntegral::model() const
{
    return (!DecayTrees_.empty() and DecayTrees_[0]) ? DecayTrees_[0]->model() : nullptr;
}

//-------------------------
const RealIntegralElement DecayTreeVectorIntegral::integral(unsigned i) const
{
    return RealIntegralElement(Diagonals_.at(i).value() * norm(DecayTrees_[i]->dataIndependentAmplitude()));
}

//-------------------------
const RealIntegralElement DecayTreeVectorIntegral::integral(unsigned i, unsigned j) const
{
    if (j < i)
        return integral(j, i);

    if (i == j)
        return integral(i);

    return RealIntegralElement(real(2. * OffDiagonals_.at(i).at(j - i - 1).value()
                                    * conj(DecayTrees_[i]->dataIndependentAmplitude())
                                    * DecayTrees_[j]->dataIndependentAmplitude()));
}

//-------------------------
const ComplexIntegralElement DecayTreeVectorIntegral::component(unsigned i, unsigned j) const
{
    if (j < i)
        return conj(component(j, i));
    if (i == j)
        return static_cast<ComplexIntegralElement>(Diagonals_.at(i));
    return OffDiagonals_.at(i).at(j - i - 1);
}

//-------------------------
DecayTreeVectorIntegral& DecayTreeVectorIntegral::operator+=(const DecayTreeVectorIntegral& rhs)
{
    if (rhs.Diagonals_.size() != Diagonals_.size())
        throw exceptions::Exception("size mismatch", "DecayTreeVectorIntegral::operator+=");

    for (size_t i = 0; i < Diagonals_.size(); ++i)
        Diagonals_[i] += rhs.Diagonals_[i];
    for (size_t i = 0; i < OffDiagonals_.size(); ++i)
        for (size_t j = 0; j < OffDiagonals_[i].size(); ++j)
            OffDiagonals_[i][j] += rhs.OffDiagonals_[i][j];
    return *this;
}

//-------------------------
DecayTreeVectorIntegral& DecayTreeVectorIntegral::operator*=(double rhs)
{
    for (size_t i = 0; i < Diagonals_.size(); ++i)
        Diagonals_[i] *= rhs;
    for (size_t i = 0; i < OffDiagonals_.size(); ++i)
        for (size_t j = 0; j < OffDiagonals_[i].size(); ++j)
            OffDiagonals_[i][j] *= rhs;
    return *this;
}

//-------------------------
DecayTreeVectorIntegral& DecayTreeVectorIntegral::reset()
{
    for (auto& elt : Diagonals_)
        elt.reset();
    for (auto& row : OffDiagonals_)
        for (auto& elt : row)
            elt.reset();
    return *this;
}

//-------------------------
const RealIntegralElementVector diagonal_integrals(const DecayTreeVectorIntegral& dtvi)
{
    RealIntegralElementVector I;
    I.reserve(dtvi.decayTrees().size());
    for (size_t i = 0; i < dtvi.decayTrees().size(); ++i)
        I.emplace_back(dtvi.integral(i));
    return I;
}

//-------------------------
const RealIntegralElementVector fit_fractions(const DecayTreeVectorIntegral& dtvi)
{
    auto I = integral(dtvi);
    auto ff = diagonal_integrals(dtvi);
    std::transform(ff.begin(), ff.end(), ff.begin(), std::bind(std::divides<RealIntegralElement>(), std::placeholders::_1, I));
    return ff;
}

//-------------------------
const RealIntegralElementVector fit_fractions(const DecayTreeVectorIntegral& dtvi, const std::vector<DecayTreeVector>& dtvv)
{
    // total integral
    auto I = integral(dtvi);

    RealIntegralElementVector ff;
    ff.reserve(dtvv.size());

    for (const auto& dtv : dtvv)
        ff.push_back(integral(dtvi, dtv) / I);
    return ff;
}

//-------------------------
const ComplexIntegralElementMatrix cached_integrals(const DecayTreeVectorIntegral& dtvi)
{
    ComplexIntegralElementMatrix I(dtvi.decayTrees().size(), ComplexIntegralElementVector(dtvi.decayTrees().size()));
    for (size_t i = 0; i < dtvi.decayTrees().size(); ++i) {
        I[i][i] = static_cast<ComplexIntegralElement>(dtvi.diagonals()[i]);
        for (size_t j = i + 1; j < dtvi.decayTrees().size(); ++j)
            I[j][i] = conj(I[i][j] = dtvi.offDiagonals()[i][j - i - 1]);
    }
    return I;
}

//-------------------------
const ComplexIntegralElementMatrix integrals(const DecayTreeVectorIntegral& dtvi)
{
    ComplexIntegralElementMatrix I(dtvi.decayTrees().size(), ComplexIntegralElementVector(dtvi.decayTrees().size()));
    for (size_t i = 0; i < dtvi.decayTrees().size(); ++i) {
        I[i][i] = static_cast<ComplexIntegralElement>(dtvi.diagonals()[i] * norm(dtvi.decayTrees()[i]->dataIndependentAmplitude()));
        for (size_t j = i + 1; j < dtvi.decayTrees().size(); ++j)
            I[j][i] = conj(I[i][j] =
                               dtvi.offDiagonals()[i][j - i - 1]
                               * conj(dtvi.decayTrees()[i]->dataIndependentAmplitude())
                               * dtvi.decayTrees()[j]->dataIndependentAmplitude());
    }
    return I;
}

//-------------------------
const RealIntegralElement integral(const DecayTreeVectorIntegral& dtvi)
{
    RealIntegralElement I(0.);
    for (unsigned i = 0; i < dtvi.diagonals().size(); ++i) {
        I += dtvi.integral(i);
        for (unsigned j = i + 1; j < dtvi.diagonals().size(); ++j)
            I += dtvi.integral(i, j);
    }
    return I;
}

//-------------------------
const RealIntegralElement integral(const DecayTreeVectorIntegral& dtvi, const DecayTreeVector& dtv)
{
    const auto& dts = dtvi.decayTrees();

    // find indices
    std::vector<unsigned> indices;
    indices.reserve(dtv.size());
    for (const auto& dt : dtv) {
        auto it = std::find(dts.begin(), dts.end(), dt);
        if (it == dts.end())
            throw exceptions::Exception("Trying to calculate integral for a DecayTree which is not in the DecayTreeVectorIntegral.", "integral");
        indices.push_back(std::distance(dts.begin(), it));
    }

    RealIntegralElement I(0.);
    for (size_t i = 0; i < indices.size(); ++i)
       for (size_t j = i; j < indices.size(); ++j)
           I += dtvi.integral(indices[i], indices[j]);

    return I;
}

}
