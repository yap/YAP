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

#ifndef yap_DecayTreeVectorIntegral_h
#define yap_DecayTreeVectorIntegral_h

#include "IntegralElement.h"

#include "fwd/DecayTreeVectorIntegral.h"

#include "fwd/DecayTree.h"
#include "fwd/Model.h"

#include <array>
#include <complex>
#include <map>

namespace yap {

/// \class DecayTreeVectorIntegral
/// \brief Stores integral components for a vector of decay trees
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
class DecayTreeVectorIntegral
{
public:

    /// constructor
    /// \param dtv DecayTreeVector to construct integral of
    explicit DecayTreeVectorIntegral(const DecayTreeVector& dtv);

    /// \return integral calculated from components
    const RealIntegralElement integral() const;

    /// \return DecayTrees_
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// \return Model this integral calculates with (via DecayTrees)
    const Model* model() const;

    /// \return Diagonals_ (const)
    const DiagonalIntegralMap& diagonals() const
    { return Diagonals_; }

    /// \return OffDiagonals_ (const)
    const OffDiagonalIntegralMap& offDiagonals() const
    { return OffDiagonals_; }

    /// grant friend status to Integrator to access components and DecayTrees_
    friend class Integrator;

private:

    /// DecayTrees to integrate
    DecayTreeVector DecayTrees_;

    /// diagonal element integrals:
    /// stores norm(dataDependentAmplitude(...)),
    /// for each DecayTree in DecayTrees_
    DiagonalIntegralMap Diagonals_;

    /// off-diagonal element integrals stores:
    /// conj([0].dataDependentAmplitude(...)) * [1].dataDependentAmplitude(...),
    /// for each pair of DecayTree's in DecayTrees_ (in the upper
    /// right triangle of the matrix of combinations)
    OffDiagonalIntegralMap OffDiagonals_;

};

/// \return vector of fit fractions of DecayTree's in DecayTreeVectorIntegral
/// \param MI DecayTreeVectorIntegral to retrieve values from
const std::vector<double> fit_fractions(const DecayTreeVectorIntegral& MI);

/// \return matrix of integral components without multiplication by free amplitudes
/// \param MI DecayTreeVectorIntegral to retrieve values from
const std::vector<std::vector<std::complex<double> > > cached_integrals(const DecayTreeVectorIntegral& MI);

/// \return matrix of integral components with multiplication by free amplitudes
/// \param MI DecayTreeVectorIntegral to retrieve values from
const std::vector<std::vector<std::complex<double> > > integrals(const DecayTreeVectorIntegral& MI);

/// \return integral of diagonal element given DiagonalMap entry
const RealIntegralElement integral(const DiagonalIntegralMap::value_type& a_A2);

/// \return integral of off-diagonal element given OffDiagonalMap entry
const RealIntegralElement integral(const OffDiagonalIntegralMap::value_type& aa_AA);

}

#endif
