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

    virtual Amp amplitude(DataPoint& d) = 0;
    virtual bool consistent() const = 0;

    unsigned int index() const {return Index_;}

protected:
    bool Recalculate_; ///

private:
    unsigned Index_; /// storage index used in DataPoint. Must be unique.

    static unsigned GlobalIndex;
};

}

#endif
