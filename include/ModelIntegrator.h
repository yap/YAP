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

#ifndef yap_ModelIntegrator_h
#define yap_ModelIntegrator_h

#include "fwd/DecayTree.h"

#include "ModelIntegral.h"

namespace yap {

/// \class ModelIntegrator
/// \brief Integrates a model
/// \author Daniel Greenwald
/// \ingroup Integration
class ModelIntegrator
{
protected:

    /// \return ModelIntegral's Diagonals_
    /// \param I ModelIntegral to access
    ModelIntegral::DiagonalMap& diagonals(ModelIntegral& I) const
    { return I.Diagonals_; }

    /// \return ModelIntegral's OffDiagonals_
    /// \param I ModelIntegral to access
    ModelIntegral::OffDiagonalMap& offDiagonals(ModelIntegral& I) const
    { return I.OffDiagonals_; }

    /// \return value stored in ModelIntegal's Diagonals_ for particular DecayTree
    /// \param I ModelIntegral to access
    /// \param dt shared_ptr to DecayTree to access diagonal integral component for
    ModelIntegral::DiagonalMap::mapped_type& diagonalComponent(ModelIntegral& I, const DecayTreeVector::value_type& dt) const
    { return I.Diagonals_.at(dt); }

    /// \return value stored in ModelIntegral's OffDiagonals_ for particular DecayTree pair
    /// \param I ModelIntegral to access
    /// \param i shared_ptr to DecayTree to access off-diagonal integral component for
    /// \param j shared_ptr to DecayTree to access off-diagonal integral component for
    ModelIntegral::OffDiagonalMap::mapped_type& offDiagonalComponent(ModelIntegral& I, const DecayTreeVector::value_type& i, const DecayTreeVector::value_type& j) const
    { return I.OffDiagonals_.at(ModelIntegral::OffDiagonalMap::key_type({i, j})); }


};

}

#endif
