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

#ifndef yap_QuantumNumbers_h
#define yap_QuantumNumbers_h

namespace yap {

/// \class QuantumNumbers
/// \brief Quantum numbers of a Particle
/// \author Johannes Rauch

class QuantumNumbers
{
public:
    QuantumNumbers(unsigned char twoJ, char P, char C, char I, char G) :
        twoJ_(twoJ), P_(P), C_(C), I_(I), G_(G) {;}
    ~QuantumNumbers() {;}

    unsigned char twoJ() const {return twoJ_;}
    double J() const {return twoJ_ * 0.5;}
    char P() const {return P_;}
    char C() const {return C_;}
    char I() const {return I_;}
    char G() const {return G_;}

    //! Checks equality of QuantumNumbers
    friend bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs);
    //! returns NOT ==
    friend bool operator!= (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
    {return !(lhs == rhs);}

private:
    unsigned char twoJ_; /// Spin * 2
    char P_; /// Parity
    char C_; /// C-parity
    char I_; /// Isospin
    char G_; /// G-parity
};

}

#endif
