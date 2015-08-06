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

#ifndef yap_Particle_h
#define yap_Particle_h

#include "DataAccessor.h"
#include "QuantumNumbers.h"

namespace yap {

/// \class Particle
/// \brief Particle base class
/// \author Johannes Rauch

/// \defgroup Particle Particle

class Particle : public DataAccessor
{
public:
    Particle();
    ~Particle();

    virtual Amp amplitude(DataPoint& d) = 0;
    virtual bool consistent() const;

    const QuantumNumbers& quantumNumbers() const {return QuantumNumbers_;}

protected:
    QuantumNumbers QuantumNumbers_;
    double Mass_;
    std::string Name_;
};

}

#endif
