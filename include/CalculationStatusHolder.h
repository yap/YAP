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

#ifndef yap_CalculationStatusHolder_h
#define yap_CalculationStatusHolder_h

#include "CalculationStatus.h"
#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

class CalculationStatusHolder
{
public:

    /// Constructor
    CalculationStatusHolder()
    {}

    /// \name calculation statuses
    /// @{

    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    virtual CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned symmetrizationIndex, unsigned dataPartitionIndex) const
    {
        DEBUG(" CalculationStatusHolder::calculationStatus " << kCalculated);
        return kCalculated;
    }

    /// @}

};

}

#endif
