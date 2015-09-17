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

#ifndef yap_SpinAmplitude_h
#define yap_SpinAmplitude_h

#include "Amp.h"
#include "AmplitudeComponentDataAccessor.h"
#include "logging.h"
#include "ParticleCombination.h"
#include "QuantumNumbers.h"

#include <array>
#include <memory>

namespace yap {

/// \class SpinAmplitude
/// \brief Abstract base class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup SpinAmplitude Spin Amplitudes

class SpinAmplitude : public AmplitudeComponentDataAccessor
{
public:

    /// Constructor
    SpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial,
                  const QuantumNumbers& final1, const QuantumNumbers& final2,
                  unsigned char twoL);

    /// Check consistency of object
    virtual bool consistent() const override;

    /// cast into string
    virtual operator std::string() const = 0;

    /// \name Getters
    /// @{

    /// Get initial QuantumNumbers
    const QuantumNumbers& initialQuantumNumbers() const
    { return InitialQuantumNumbers_; }

    /// Get QuantumNumbers of daughters const
    const std::array<QuantumNumbers, 2>& finalQuantumNumbers() const
    { return FinalQuantumNumbers_; }

    /// Get relative angular momentum between daughters * 2
    unsigned char twoL() const
    { return TwoL_; }

    /// Get relative angular momentum between daughters
    double L() const
    { return 0.5 * TwoL_; }

    /// @}

    /// Compare SpinAmplitude objects
    friend bool operator== (const SpinAmplitude& lhs, const SpinAmplitude& rhs)
    { return typeid(lhs) == typeid(rhs) && lhs.equals(rhs); }

protected:

    /// \return Complex spin amplitude evaluated at data point
    /// \param d DataPoint to evaluate on
    virtual Amp calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) override = 0;

    /// Check if SpinAmplitudes are equal
    virtual bool equals(const SpinAmplitude& rhs) const;

    /// Initial-state quantum numbers
    QuantumNumbers InitialQuantumNumbers_;

    /// array of final-state quantum numbers
    std::array<QuantumNumbers, 2> FinalQuantumNumbers_;

    /// relative angular momentum between daughters * 2
    unsigned char TwoL_;

};

struct SharedSpinAmplitudeComparator {
/// Compare SpinAmplitude shared_ptr's
    bool operator() (const std::shared_ptr<SpinAmplitude>& lhs, const std::shared_ptr<SpinAmplitude>& rhs) const
    { return lhs.get() == rhs.get() || *lhs == *rhs; }

};


}

#endif
