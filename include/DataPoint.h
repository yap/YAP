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

#ifndef yap_DataPoint_h
#define yap_DataPoint_h

#include <Constants.h>
#include <DataAccessor.h>
#include <FourVector.h>

#include <complex>
#include <memory>
#include <vector>

namespace yap {

class CachedDataValue;
class FourMomenta;
class HelicityAngles;
class MeasuredBreakupMomenta;
class ParticleCombination;

/// \class DataPoint
/// \brief Class for holding data and cached values per data point for fast calculation
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Data Data-related classes

class DataPoint
{
public:

    /// \name Constructors
    /// @{

    DataPoint(const DataAccessorSet& S);

    /// @}

    /// \return number of data accessor rows
    size_t nDataAccessors() const
    { return Data_.size(); }

    /// \return number of sym indices rows for data accessor
    /// \param i index of DataAccessor
    size_t nSymIndices(unsigned i) const
    { return Data_[i].size(); }

    /// \return number of elements for data accessor
    /// \param i index of DataAccessor
    /// \param j index of symmetrization
    size_t nElements(unsigned i, unsigned j = 0) const
    { return Data_[i][j].size(); }

    /// \return size of data point
    unsigned dataSize() const;

    /// \return string of size of data point
    std::string dataSizeString() const
    { return "Size of DataPoint: " + std::to_string(dataSize()) + " byte (for " + std::to_string(Data_.size()) + " data accessors"; }

    /// check that two DataPoint's have same internal structure
    friend bool equalStructure(const DataPoint& A, const DataPoint& B);

    /// check that two DataPoint's are equal
    friend bool operator==(const DataPoint& lhs, const DataPoint& rhs)
    { return lhs.Data_ == rhs.Data_; }

    /// grant friend status to CachedDataValue to access Data_
    friend class CachedDataValue;

private:

    /// Data storage for all DataAccessors
    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<std::vector<double> > > Data_;

};

}

#endif
