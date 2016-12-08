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

#include <vector>

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
    static std::vector<ModelComponentIntegral>& integrals(ModelIntegral& I);

    /// \return DecayTreeVectorIntegral's Diagonals_
    /// \param I DecayTreeVectorIntegral to access
    static RealIntegralElementVector& diagonals(DecayTreeVectorIntegral& I);

    /// \return DecayTreeVectorIntegral's OffDiagonals_
    /// \param I DecayTreeVectorIntegral to access
    static ComplexIntegralElementMatrix& offDiagonals(DecayTreeVectorIntegral& I);

    /// zero-out a DecayTreeVectorIntegral
    static DecayTreeVectorIntegral& reset(DecayTreeVectorIntegral& I);

};

}

#endif
