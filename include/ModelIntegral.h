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

#include "fwd/ModelIntegral.h"

#include "fwd/DecayTree.h"
#include "fwd/IntegralElement.h"
#include "fwd/Model.h"

#include "DecayTreeVectorIntegral.h"

namespace yap {

/// \class ModelIntegral
/// \brief Stores integral of a full model
/// \author Daniel Greenwald
/// \ingroup Integration
class ModelIntegral
{

public:

    /// constructor
    /// \param model Model to integrate
    ModelIntegral(const Model& model);

    const IntegralMap& integrals() const
    { return Integrals_; }

    /// grant friend status to Integrator for nonconst access Integrals_
    friend class Integrator;

private:

    /// Map of (admixture factor -> component integral)
    IntegralMap Integrals_;

};

/// \return integral calculated from components
const RealIntegralElement integral(const ModelIntegral& MI);


}

#endif
