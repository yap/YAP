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

#include "fwd/DataPartition.h"
#include "fwd/DecayTreeVectorIntegral.h"

#include "DecayTreeVectorIntegrator.h"

namespace yap {

/// \class ImportanceSampler
/// \brief Calculates DecayTreeVectorIntegral using importance sampling
/// \author Daniel Greenwald
/// \ingroup Integration
class ImportanceSampler : public DecayTreeVectorIntegrator
{

public:

    /// Update calculation of DecayTreeVectorIntegral
    static void calculate(DecayTreeVectorIntegral& I, DataPartitionVector& DPV);

    /// Update calculation of DecayTreeVectorIntegral
    static void calculate(DecayTreeVectorIntegral& I, DataPartition& D);

private:

    /// perform partial calculation for one data partition
    static unsigned partialCalculation(DecayTreeVectorIntegral& I, DataPartition& D);

};

}

#endif
