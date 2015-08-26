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

#ifndef yap_CanonicalSpinAmplitude_h
#define yap_CanonicalSpinAmplitude_h

#include "SpinAmplitude.h"

namespace yap {

/// \class CanonicalSpinAmplitude
/// \brief Class implementing a canonical spin amplitude, i.e. with defined relative angular momentum.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude

class CanonicalSpinAmplitude : public SpinAmplitude
{
public:

    /// Constructor
    CanonicalSpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char l);

    /// \return Complex spin amplitude evaluated at data point
    /// \param d DataPoint to evaluate on
    virtual Amp amplitude(DataPoint& d) override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// cast into string
    operator std::string() const override;

    /// \name Getters
    /// @{

    /// Get relative angular momentum between daughters
    unsigned char decayAngularMomentum() const
    { return L_; }

    /// @}


private:

    /// Check if SpinAmplitudes are equal
    bool equals(const SpinAmplitude& rhs) const override;

    /// relative angular momentum between daughters
    unsigned char L_;

};

}

#endif
