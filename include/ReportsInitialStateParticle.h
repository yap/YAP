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

#ifndef yap_ReportsInitialStateParticle_h
#define yap_ReportsInitialStateParticle_h

namespace yap {

class InitialStateParticle;

/// \name ReportsInitialStateParticle
/// \brief Base class for all classes that report which InitialStateParticle they belong
/// \author Johannes Rauch, Daniel Greenwald

class ReportsInitialStateParticle
{
public:

    /// get raw pointer to initial state particle
    virtual InitialStateParticle* initialStateParticle() = 0;

    /// get raw pointer to initial state particle (const)
    const InitialStateParticle* initialStateParticle() const
    { return const_cast<ReportsInitialStateParticle*>(this)->initialStateParticle(); }

    /// Default constructor
    ReportsInitialStateParticle() {}

    /// virtual destructor
    virtual ~ReportsInitialStateParticle() = default;

    /// default copy constructor
    ReportsInitialStateParticle(const ReportsInitialStateParticle& other) = default;

    /// default move constructor
    ReportsInitialStateParticle(ReportsInitialStateParticle&& other) = default;

    /// default copy assignment operator
    ReportsInitialStateParticle& operator=(const ReportsInitialStateParticle& rhs) = default;

    /// default move assignment operator
    ReportsInitialStateParticle& operator=(ReportsInitialStateParticle&& rhs) = default;

};

}

#endif
