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

#include <array>
#include <complex>
#include <map>

namespace yap {

class ModelIntegral
{
public:

    /// constructor
    /// \param dtv DecayTreeVector to construct integral of
    ModelIntegral(const DecayTreeVector& dtv);

    const double integral() const;

    /// \name typedefs
    /// @{

    /// \typedef map type for diagonal integral elements
    using DiagonalMap = std::map<DecayTreeVector::value_type, double>;

    /// \typdef map type for off-diagonal integral elements
    using OffDiagonalMap = std::map<std::array<DecayTreeVector::value_type, 2>, std::complex<double> >;

    /// @}

private:

    /// vector of DecayTree's to be used to calculate integral
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

/// \return integral of diagonal element given DiagonalMap entry
const double integral(const ModelIntegral::DiagonalMap::value_type& a_A2);

/// \return integral of off-diagonal element given OffDiagonalMap entry
const double integral(const ModelIntegral::OffDiagonalMap::value_type& aa_AA);

}

#endif
