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

#include <ostream>

namespace yap {

/// \class QuantumNumbers
/// \brief Quantum numbers of a Particle
/// \author Johannes Rauch

class QuantumNumbers
{
public:
    /// Constructor
    QuantumNumbers(unsigned char twoJ, char P, char C, char I, char G);

    /// \name Getters
    /// @{

    /// return spin * 2
    unsigned char twoJ() const {return twoJ_;}

    /// return spin
    double J() const {return twoJ_ * 0.5;}

    /// return parity
    char P() const {return P_;}

    /// return C-parity
    char C() const {return C_;}

    /// return isospin
    char I() const {return I_;}

    /// return G-parity
    char G() const {return G_;}

    /// @}

    //! Check equality of QuantumNumbers
    friend bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs);

    //! returns NOT ==
    friend bool operator!= (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
    { return !(lhs == rhs); }

private:
    unsigned char twoJ_; /// Spin * 2
    char P_; /// Parity
    char C_; /// C-parity
    char I_; /// Isospin
    char G_; /// G-parity
};

/// Overload << operator
std::ostream& operator<< (std::ostream&, const QuantumNumbers&);

}

#endif
