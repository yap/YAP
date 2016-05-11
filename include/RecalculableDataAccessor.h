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

#ifndef yap_RecalculableDataAccessor_h
#define yap_RecalculableDataAccessor_h

#include "fwd/RecalculableDataAccessor.h"
#include "fwd/DataPartition.h"

#include "DataAccessor.h"

namespace yap {

class RecalculableDataAccessor : public DataAccessor
{
public:
    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    RecalculableDataAccessor(ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer)
        : DataAccessor(equiv) {}

    /// calculate for every data point in a data partition
    virtual void calculate(DataPartition& D) const = 0;

};

}

#endif
