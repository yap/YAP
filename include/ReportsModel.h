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

#ifndef yap_ReportsModel_h
#define yap_ReportsModel_h

namespace yap {

class Model;

/// \name ReportsModel
/// \brief Base class for all classes that report which Model they belong
/// \author Johannes Rauch, Daniel Greenwald

class ReportsModel
{
public:

    /// get raw pointer to Model (const)
    virtual const Model* model() const = 0;

    /// get raw pointer to Model
    Model* model()
    { return const_cast<Model*>(const_cast<const ReportsModel*>(this)->model()); }

    /// Default constructor
    ReportsModel() {}

    /// virtual destructor
    virtual ~ReportsModel() = default;

    /// default copy constructor
    ReportsModel(const ReportsModel& other) = default;

    /// default move constructor
    ReportsModel(ReportsModel&& other) = default;

    /// default copy assignment operator
    ReportsModel& operator=(const ReportsModel& rhs) = default;

    /// default move assignment operator
    ReportsModel& operator=(ReportsModel&& rhs) = default;

};

}

#endif
