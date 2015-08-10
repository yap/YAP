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

#ifndef yap_SpinAmplitude_h
#define yap_SpinAmplitude_h

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "QuantumNumbers.h"

#include <array>

namespace yap {

/// \class InitialStateParticle
/// \brief Class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald

class SpinAmplitude : public AmplitudeComponent, public DataAccessor
{
public:
    /// Constructor
    SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2);

    virtual Amp amplitude(DataPoint& d) override;
    virtual bool consistent() const override;

    /// \name Getters
    /// @{

    /// Get initial QuantumNumbers
    const QuantumNumbers& initialQuantumNumbers() const {return InitialQuantumNumbers_;}

    /// Get final QuantumNumbers of 1st daughter
    const QuantumNumbers& finalQuantumNumbersA() const {return FinalQuantumNumbers_[0];}

    /// Get initial QuantumNumbers of 2nd daughter
    const QuantumNumbers& finalQuantumNumbersB() const {return FinalQuantumNumbers_[1];}

    /// @}

private:
    QuantumNumbers InitialQuantumNumbers_;
    std::array<QuantumNumbers, 2> FinalQuantumNumbers_;

};

}

#endif
