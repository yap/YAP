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

#ifndef yap_ImportanceSampler_h
#define yap_ImportanceSampler_h

#include "Integrator.h"

#include "fwd/DataPartition.h"
#include "fwd/DecayTreeVectorIntegral.h"
#include "fwd/ModelIntegral.h"

namespace yap {

/// \class ImportanceSampler
/// \brief Calculates DecayTreeVectorIntegral using importance sampling
/// \author Daniel Greenwald
/// \ingroup Integration
class ImportanceSampler : public Integrator
{

public:

    /// Update calculation of ModelIntegral
    static void calculate(ModelIntegral& I, DataPartitionVector& DPV);

    /// Update calculation of ModelIntegral
    static void calculate(ModelIntegral& I, DataPartition& D);

private:

    /// \return integral_sub_map for all changed trees
    static std::vector<DecayTreeVectorIntegral*> select_changed(ModelIntegral& I);

    /// perform partial calculation of one integral component for one data partition
    static unsigned partially_calculate_subIntegral(DecayTreeVectorIntegral& I, DataPartition& D);

    /// perform partial calculation for one data partition
    static unsigned partially_calculate(std::vector<DecayTreeVectorIntegral*>& J, DataPartition& D);

};

}

#endif
