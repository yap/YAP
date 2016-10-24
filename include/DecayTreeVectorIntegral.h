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
///   o Diagonal components: one for each DecayTree, stored as |A|^2
///     for the tree and returned as |a|^2 * stored value when
///     integral(i) is called
///   o Off-diagonal components: one for each pair of DecayTrees (i, j),
///     stored as conj(A_i) * A_j and returned as
///     2 * real(conj(a_i) * a_i * stored value when integral(i, j)
///     is called.
class DecayTreeVectorIntegral
{
public:

    /// constructor
    /// \param dtv DecayTreeVector to construct integral of
    explicit DecayTreeVectorIntegral(const DecayTreeVector& dtv);

    /// \return integral of diagonal component
    const RealIntegralElement integral(unsigned i) const;

    /// \return integral of off-diagonal components := (i,j) + (j,i)*, if i != j; else (i,i)
    const RealIntegralElement integral(unsigned i, unsigned j) const;

    /// \return DecayTrees_
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// \return Model this integral calculates with (via DecayTrees)
    const Model* model() const;

    /// \return Diagonals_ (const)
    const RealIntegralElementVector& diagonals() const
    { return Diagonals_; }

    /// \return OffDiagonals_ (const)
    const ComplexIntegralElementMatrix& offDiagonals() const
    { return OffDiagonals_; }

    /// casts diagonal components into off-diagonal type, conjugates
    /// upper-triangle off-diagonals to return lower-triangle members.
    /// \return (copy of) component---integral of conj(A_i) * A_j
    const ComplexIntegralElement component(unsigned i, unsigned j) const;

    // addition operator
    DecayTreeVectorIntegral& operator+=(const DecayTreeVectorIntegral& rhs);

    // multiplication operator
    DecayTreeVectorIntegral& operator*=(double rhs);

    /// grant friend status to Integrator to access components and DecayTrees_
    friend class Integrator;

protected:

    DecayTreeVectorIntegral& reset();

private:

    /// DecayTrees to integrate
    DecayTreeVector DecayTrees_;

    /// diagonal element integrals:
    /// stores norm(dataDependentAmplitude(...)),
    /// for each DecayTree in DecayTrees_
    RealIntegralElementVector Diagonals_;

    /// off-diagonal element integrals stores:
    /// conj([0].dataDependentAmplitude(...)) * [1].dataDependentAmplitude(...),
    /// for each pair of DecayTree's in DecayTrees_ (in the upper
    /// right triangle of the matrix of combinations)
    ComplexIntegralElementMatrix OffDiagonals_;

};

/// \return integral calculated from components of DecayTreeVectorIntegral
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
const RealIntegralElement integral(const DecayTreeVectorIntegral& dtvi);

/// \return integral calculated from given components of DecayTreeVectorIntegral
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
/// \param dtv DecayTrees to integrate over
const RealIntegralElement integral(const DecayTreeVectorIntegral& dtvi, const DecayTreeVector& dtv);

/// \return matrix of integral components without multiplication by free amplitudes
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
const ComplexIntegralElementMatrix cached_integrals(const DecayTreeVectorIntegral& dtvi);

/// \return matrix of integral components with multiplication by free amplitudes
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
const ComplexIntegralElementMatrix integrals(const DecayTreeVectorIntegral& dtvi);

/// \return diagonal integrals
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
const RealIntegralElementVector diagonal_integrals(const DecayTreeVectorIntegral& dtvi);

/// \return vector of fit fractions of DecayTree's in DecayTreeVectorIntegral
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
const RealIntegralElementVector fit_fractions(const DecayTreeVectorIntegral& dtvi);

/// \return vector of fit fractions of groups of DecayTree's in DecayTreeVectorIntegral
/// \param dtvi DecayTreeVectorIntegral to retrieve values from
/// \param dtvv a fit fraction will be calculated for each DecayTreeVector
const RealIntegralElementVector fit_fractions(const DecayTreeVectorIntegral& dtvi, const std::vector<DecayTreeVector>& dtvv);

}

#endif
