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

#ifndef yap__ParameterSet_h
#define yap__ParameterSet_h

#include <initializer_list>
#include <vector>

namespace yap {

/// \class ParameterSet
/// \brief Base class for all objects holding fitable parameters
/// \author Johannes Rauch, Daniel Greenwald
///
/// Manages locking and unlocking of parameters

class ParameterSet
{
public:

    /// \enum ParameterStatus
    enum ParameterStatus {
        kChanged   = -1,        ///< Parameter is free and has been changed
        kFixed     = 0,         ///< Parameter is fixed
        kUnchanged = +1,        ///< Parameter is free but has not been changed
    };

    /// \name Constructors and operators
    /// @{

    /// constructor
    /// \param pars An initializer_list of doubles to set into the parameters
    /// \param status the status to set for all parameters
    ParameterSet(std::initializer_list<double> pars = {}, ParameterStatus status = kFixed)
        : Parameters_(pars), ParameterStatuses_(pars.size(), status)
    {}

    /// assignment operator
    ParameterSet& operator=(std::initializer_list<double> pars);

    /// @}

    /// Get parameters
    /*std::vector<double>& parameters()
    { return Parameters_; }*/

    /// Get (const) parameters
    const std::vector<double>& parameters() const
    { return Parameters_; }

    /// Get parameter statuses
    /*std::vector<ParameterStatus>& parameterStatuses()
    { return ParameterStatuses_; }*/

    /// Get (const) parameters
    const std::vector<ParameterStatus>& parameterStatuses() const
    { return ParameterStatuses_; }

    /// Check consistency of object
    virtual bool consistent() const
    { return Parameters_.size() == ParameterStatuses_.size(); }

protected:

    /// synchronize #ParameterStatuses_ to #Parameters_
    /// \param status Status to initialize new parameters with (default = #kFixed)
    void synchronizeParameterStatuses(ParameterStatus status = kFixed)
    { ParameterStatuses_.resize(Parameters_.size(), status); }

    /// Parameters
    std::vector<double> Parameters_;

    /// Status of parameters
    std::vector<ParameterStatus> ParameterStatuses_;

};

}

#endif
