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

#ifndef yap_DecayTree_h
#define yap_DecayTree_h

#include <CachedDataValue.h>
#include <Parameter.h>

#include <vector>

namespace yap {

class CachedDataValue;
class Parameter;

/// \class DecayTree
/// \brief Class holding vectors of fixed and free amplitudes that define a decay tree
/// \author Johannes Rauch, Daniel Greenwald
class DecayTree
{
public:

    /// default constructor
    DecayTree() {}

    /// constructor
    DecayTree(const CachedDataValueVector& fixedAmplitudes, const ComplexParameterVector& freeAmplitudes) :
        fixedAmplitudes_(fixedAmplitudes), freeAmplitudes_(freeAmplitudes)
    {}

    /// \name Getters
    /// @{

    const CachedDataValueVector& fixedAmplitudes() const
    { return fixedAmplitudes_; }

    CachedDataValueVector fixedAmplitudes()
    { return fixedAmplitudes_; }

    const ComplexParameterVector& freeAmplitudes() const
    { return freeAmplitudes_; }

    ComplexParameterVector freeAmplitudes()
    { return freeAmplitudes_; }

    /// @}

private:

    CachedDataValueVector fixedAmplitudes_;
    ComplexParameterVector freeAmplitudes_;

};

}

#endif
