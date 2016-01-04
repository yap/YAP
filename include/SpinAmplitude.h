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

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "QuantumNumbers.h"
#include "DecayChannel.h"

#include <array>
#include <cstdlib>
#include <memory>

namespace yap {

/// \class SpinAmplitude
/// \brief Abstract base class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup SpinAmplitude Spin Amplitudes
class SpinAmplitude : public AmplitudeComponent, public DataAccessor
{
public:

    /// \return Whether angular momentum is conserved in J -> j1 + j2 with orbital angular momentum l
    /// \param two_J 2 * spin of initial state
    /// \param two_j1 2 * spin of first daughter
    /// \param two_j2 2 * spin of second daughter
    /// \param l orbital angular momentum
    static constexpr bool conserves(unsigned two_J, unsigned two_j1, unsigned two_j2, int l)
    {
        // check that the spins are consistent, and that the triangle requirements for (Jlj) and (j1j2j) can be met simultaneously
        return is_even(two_J + two_j1 + two_j2)
               and (std::min(two_j1 + two_j2, two_J + 2 * l) >= std::max(std::abs<int>(two_j1 - two_j2), std::abs<int>(two_J - 2 * l)));
    }

    /// Constructor
    /// \param intial quantum numbers of Initial-state
    /// \param final1 quantum numbers of first daughter
    /// \param final2 quantum numbers of second daughter
    /// \param orbital angular momentum
    SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned l);

    /// Check consistency of object
    virtual bool consistent() const override;

    /// cast into string
    virtual operator std::string() const;

    /// \name Getters
    /// @{

    /// Get initial QuantumNumbers
    const QuantumNumbers& initialQuantumNumbers() const
    { return InitialQuantumNumbers_; }

    /// Get QuantumNumbers of daughters const
    const std::array<QuantumNumbers, 2>& finalQuantumNumbers() const
    { return FinalQuantumNumbers_; }

    /// Get orbital angular momentum
    unsigned l() const
    { return L_; }

    /// include const access to ISP
    using BelongsToInitialStateParticle::initialStateParticle;

    /// Get raw pointer to InitialStateParticle
    InitialStateParticle* initialStateParticle() override
    { return InitialStateParticle_; }

    /// @}

    /// Add symmetrization indices for ParticleCombination.
    /// Must be overloaded in derived classes to both conditionally add ParticleCombination's
    /// and to alter and multiply them (to accomodate helicity, for example)
    virtual ParticleCombinationVector addSymmetrizationIndices(std::shared_ptr<const ParticleCombination> pc) = 0;

    /// grant friend access to DecayChannel to set InitialStateParticle
    friend DecayChannel;

protected:

    /// set raw pointer to owning InitialStateParticle
    virtual void setInitialStateParticle(InitialStateParticle* isp)
    { InitialStateParticle_ = isp; }

    /// check equality
    virtual bool equals(const SpinAmplitude& other) const;

private:

    /// raw pointer to owning InitialStateParticle
    InitialStateParticle* InitialStateParticle_;

    /// Initial-state quantum numbers
    QuantumNumbers InitialQuantumNumbers_;

    /// array of final-state quantum numbers
    std::array<QuantumNumbers, 2> FinalQuantumNumbers_;

    /// orbital angular momentum
    unsigned L_;

    /// equality operator
    friend bool operator==(const SpinAmplitude& A, const SpinAmplitude& B)
    { return typeid(A) == typeid(B) and A.equals(B); }

};

/// convert to string
inline std::string to_string(const SpinAmplitude& sa)
{ return (std::string)sa; }

/// << operator
inline std::ostream& operator<< (std::ostream& os, const SpinAmplitude& sa)
{ os << to_string(sa); return os; }

struct SharedSpinAmplitudeComparator {
/// Compare SpinAmplitude shared_ptr's
    bool operator() (const std::shared_ptr<SpinAmplitude>& lhs, const std::shared_ptr<SpinAmplitude>& rhs) const
    { return lhs.get() == rhs.get() or * lhs == *rhs; }

};


}

#endif
