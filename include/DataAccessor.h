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

#ifndef yap_DataAccessor_h
#define yap_DataAccessor_h

#include "Amp.h"
#include "DataPoint.h"
#include "ParticleCombination.h"

#include <map>

namespace yap {

/// \name DataAccessor
/// \brief Base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald

class DataAccessor
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Default constructor
    DataAccessor();

    /// Copy constructor
    DataAccessor(const DataAccessor& other);

    // Defaulted move constructor
    // Defaulted destructor
    // Defaulted move assignment operator

    /// @}

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    unsigned index() const
    { return Index_; }

    /// \return index inside row of DataPoint for the requested symmetrization
    unsigned symmetrizationIndex(std::shared_ptr<ParticleCombination> c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return list of all ParticleCombinations
    std::vector<std::shared_ptr<ParticleCombination> > particleCombinations() const;

    /// @}

    bool consistent() const;

    /// \name Symmetrization functions
    /// @{

    /// add symmetrizationIndex to SymmetrizationIndices_
    void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c);

    /// checks equality of symmetrizations, particular to the
    /// DataAccessor object. Can be overloaded to reduce redundant calculations
    virtual bool areEqual(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B)
    { return A == B; }

    /// @}

protected:

    /// Flag to mark whether recalculation needs to take place
    bool Recalculate_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    std::map<std::shared_ptr<ParticleCombination>, unsigned, std::owner_less<std::shared_ptr<ParticleCombination> > > SymmetrizationIndices_;

private:

    /// storage index used in DataPoint. Must be unique.
    unsigned Index_;

    /// static counter for setting indices
    static unsigned GlobalIndex;
};

}

#endif
