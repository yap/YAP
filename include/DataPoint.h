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

    /// 4-momenta constructor
    /// \param P vector of four-momenta to initialize to
    DataPoint(const std::vector<FourVector<double> >& P) : FSPFourMomenta_(P) {}

    /// initializes fsp four momenta to #FourVector_0
    /// \param n number of final state particles
    DataPoint(unsigned n) : FSPFourMomenta_(n, FourVector_0) {}

    // /// Invariant mass constructor
    //DataPoint(const std::map<std::shared_ptr<ParticleCombination>, double>& m2);

    /// @}

    /// get fsp four momenta (const)
    const std::vector<FourVector<double> >& finalStateFourMomenta() const
    { return FSPFourMomenta_; }

    /// get non-fsp four momenta (const)
    const std::vector<FourVector<double> >& fourMomenta() const
    { return FourMomenta_; }

    /// set four momenta into data point
    /// \param fourMomenta Four-momenta to set to
    /// \param check Whether to check equality of old and new four-momenta
    /// \return Whether new momenta == old momenta if check == true, else false
    bool setFinalStateFourMomenta(const std::vector<FourVector<double> >& fourMomenta, bool check = true);

    /// print information about the size of the DataPoint and its members
    void printDataSize();

    /// reserve space in Data_
    void allocateStorage(std::shared_ptr<FourMomenta> fourMom, const DataAccessorSet& dataAccessors);

    friend class CachedDataValue;
    friend class DataAccessor;
    friend class DataSet;
    friend class FourMomenta;

private:

    /// vector of 4-momenta of final-state particles in event
    std::vector<FourVector<double> > FSPFourMomenta_;

    /// Vector of 4-momenta of non-final-state particles in event
    std::vector<FourVector<double> > FourMomenta_;

    /// Data storage for all DataAccessors
    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<std::vector<double> > > Data_;

};

}

#endif
