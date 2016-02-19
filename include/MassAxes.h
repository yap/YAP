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

#ifndef yap_MassAxes_h
#define yap_MassAxes_h

#include "ParticleCombination.h"

namespace yap {

/// \class MassAxes
/// \brief ParticleCombinationVector specialized to contain axes for defining a phase-space coordinate
/// \author Daniel Greenwald
class MassAxes : public ParticleCombinationVector
{
public:

    /// Default constructor
    MassAxes() : ParticleCombinationVector() {}

    /// grant friend status to Model to create MassAxes
    friend class Model;

protected:

    /// protected constructor, one must use FourMomenta::getMassAxes
    MassAxes(const ParticleCombinationVector& axes) : ParticleCombinationVector(axes) {}

};

/// convert to string
inline std::string to_string(const MassAxes& A)
{
    std::string s = "";
    for (auto a : A)
        s += indices_string(*a) + ", ";
    s.erase(s.size() - 2, 2);
    return s;
}

}

#endif
