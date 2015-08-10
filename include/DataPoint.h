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

namespace yap {

/// \class DataPoint
/// \brief Class for holding data and cached values per data point for fast calculation
/// \author Johannes Rauch, Daniel Greenwald

class DataPoint
{
public:

    /// Default constructor
    DataPoint();

    // /// 4-momenta constructor
    // DataPoint(const std::vector<TLorentzVector>& P);

    // /// Invariant mass constructor
    // DataPoint(const std::vector<double>& S);

protected:

    // /// Actual data values
    // /// first index is for the DataAccessor
    // /// second index is for the symmeterization state (as known by the DataAccessor)
    // /// third index is internal to the DataAccessor
    // std::vector<std::vector<std::vector<double> > > Data_;

    // std::vector<TLorentzVector> FourMomentum_;
    // std::vector<double> SquaredMass_;

};

}

#endif
