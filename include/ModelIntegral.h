/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_ModelIntegral_h
#define yap_ModelIntegral_h

#include "fwd/DecayTree.h"
#include "fwd/Model.h"
#include "fwd/ModelIntegral.h"

#include <array>
#include <complex>
#include <map>

namespace yap {

/// \class IntegralElement
/// \brief Holds the values of a component of an integral
/// \author Daniel Greenwald
/// \ingroup Integration
template <typename T>
struct IntegralElement {
    /// integral value
    T value;

    /// constructor
    /// \param val initial value of integral component
    IntegralElement(T val = 0) : value(val) {}
};

/// \class ModelIntegral
/// \brief Stores integral components for a model
/// \author Daniel Greenwald
/// \defgroup Integration Classes related to model integration
///
/// Each DecayTree holds a product of free amplitudes ("a", below)
/// and a product of data-dependent amplitudes ("A", below)
///
/// This class holds two types of components:\n
///   o Diagonal components: one for each DecayTree,
///     stored as |A|^2 for the tree and returned as |a|^2 * stored value
///   o Off-diagonal components: one for each pair of DecayTrees (i, j),
///     stored as conj(A_i) * A_j and returned as 2 * real(conj(a_i) * a_i * stored value)
class ModelIntegral
{
public:

    /// constructor
    /// \param dtv DecayTreeVector to construct integral of
    ModelIntegral(const DecayTreeVector& dtv);

    /// \return integral calculated from components
    const RealIntegralElement integral() const;

    /// \return DecayTrees_
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// \return Model this integral calculates with (via DecayTrees)
    const Model* model() const;

    /// \name typedefs
    /// @{

    /// \typedef map type for diagonal integral elements
    using DiagonalMap = std::map<DecayTreeVector::value_type, RealIntegralElement>;

    /// \typdef map type for off-diagonal integral elements
    using OffDiagonalMap = std::map<std::array<DecayTreeVector::value_type, 2>, ComplexIntegralElement>;

    /// @}

    /// \return Diagonals_ (const)
    const ModelIntegral::DiagonalMap& diagonals() const
    { return Diagonals_; }

    /// \return OffDiagonals_ (const)
    const ModelIntegral::OffDiagonalMap& offDiagonals() const
    { return OffDiagonals_; }

    /// grant friend status to ModelIntegrator to access components and DecayTrees_
    friend class ModelIntegrator;

private:

    /// DecayTrees to integrate
    DecayTreeVector DecayTrees_;

    /// diagonal element integrals:
    /// stores norm(dataDependentAmplitude(...)),
    /// for each DecayTree in DecayTrees_
    DiagonalMap Diagonals_;

    /// off-diagonal element integrals stores:
    /// conj([0].dataDependentAmplitude(...)) * [1].dataDependentAmplitude(...),
    /// for each pair of DecayTree's in DecayTrees_ (in the upper
    /// right triangle of the matrix of combinations)
    OffDiagonalMap OffDiagonals_;

};

/// \return addition of two RealIntegralElements
inline const RealIntegralElement operator+(const RealIntegralElement& A, const RealIntegralElement& B)
{ return RealIntegralElement(A.value + B.value); }

/// \return integral of diagonal element given DiagonalMap entry
const RealIntegralElement integral(const ModelIntegral::DiagonalMap::value_type& a_A2);

/// \return integral of off-diagonal element given OffDiagonalMap entry
const RealIntegralElement integral(const ModelIntegral::OffDiagonalMap::value_type& aa_AA);

}

#endif
