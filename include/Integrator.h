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

#ifndef yap_Integrator_h
#define yap_Integrator_h

#include "fwd/DecayTree.h"
#include "fwd/DecayTreeVectorIntegral.h"
#include "fwd/ModelIntegral.h"

namespace yap {

/// \class Integrator
/// \brief Integrates a model
/// \author Daniel Greenwald
/// \ingroup Integration
class Integrator
{
protected:

    /// \return ModelIntegral's IntegralMap
    /// \param I ModelIntegral to access
    static IntegralMap& integrals(ModelIntegral& I);

    /// \return DecayTreeVectorIntegral's Diagonals_
    /// \param I DecayTreeVectorIntegral to access
    static DiagonalIntegralMap& diagonals(DecayTreeVectorIntegral& I);

    /// \return DecayTreeVectorIntegral's OffDiagonals_
    /// \param I DecayTreeVectorIntegral to access
    static OffDiagonalIntegralMap& offDiagonals(DecayTreeVectorIntegral& I);

    /// \return value stored in ModelIntegal's Diagonals_ for particular DecayTree
    /// \param I DecayTreeVectorIntegral to access
    /// \param dt shared_ptr to DecayTree to access diagonal integral component for
    static DiagonalIntegralMap::mapped_type& diagonalComponent(DecayTreeVectorIntegral& I, const DecayTreeVector::value_type& dt);

    /// \return value stored in DecayTreeVectorIntegral's OffDiagonals_ for particular DecayTree pair
    /// \param I DecayTreeVectorIntegral to access
    /// \param i shared_ptr to DecayTree to access off-diagonal integral component for
    /// \param j shared_ptr to DecayTree to access off-diagonal integral component for
    static OffDiagonalIntegralMap::mapped_type& offDiagonalComponent(DecayTreeVectorIntegral& I, const DecayTreeVector::value_type& i, const DecayTreeVector::value_type& j);

};

}

#endif