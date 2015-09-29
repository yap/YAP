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

#ifndef yap_QuantumNumbers_h
#define yap_QuantumNumbers_h

#include <cmath>
#include <ostream>

namespace yap {

/// \class QuantumNumbers
/// \brief Quantum numbers of a Particle
/// \author Johannes Rauch

class QuantumNumbers
{
public:

    /// \name Constructors
    /// @{

    /// Constructor
    QuantumNumbers(unsigned char twoJ, char P, char C, unsigned char twoI, char G, char Q);

    /// IJPQ(lambda) constructor
    QuantumNumbers(unsigned char twoI, unsigned char twoJ, char P, char Q)
        : QuantumNumbers(twoJ, P, 0, twoI, 0, Q) {}

    /// Default constructor
    /// is inconsistent
    QuantumNumbers()
        : QuantumNumbers(0, 0, 0, 0) {}

    /// @}

    /// check consistency
    virtual bool consistent() const;

    /// \name Getters
    /// @{

    /// \return spin * 2
    unsigned char twoJ() const
    { return TwoJ_; }

    /// \return spin
    double J() const
    { return TwoJ_ * 0.5; }

    /// \return parity
    signed char P() const
    { return P_; }

    /// \return C-parity
    signed char C() const
    { return C_; }

    /// \return Isospin * 2
    unsigned char twoI() const
    { return TwoI_; }

    /// \return Isospin
    double I() const
    { return TwoI_ * 0.5; }

    /// \return G-parity
    signed char G() const
    { return G_; }

    /// \return Electric charge
    char Q() const
    { return Q_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set Spin
    void setJ(double J)
    { TwoJ_ = std::round(2.*J); }

    /// Set 2 * Spin
    void setTwoJ(unsigned char J)
    { TwoJ_ = J; }

    /// @}

    /// equality operator
    friend bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs);

    /// returns NOT ==
    friend bool operator!= (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
    { return !(lhs == rhs); }

private:

    /// Spin * 2
    unsigned char TwoJ_;

    /// Parity
    signed char P_;

    /// C-parity
    signed char C_;

    /// Isospin * 2
    unsigned char TwoI_;

    /// G-parity
    signed char G_;

    /// Electric charge
    char Q_;

};

/// Overload << operator
std::ostream& operator<< (std::ostream&, const QuantumNumbers&);

}

#endif
