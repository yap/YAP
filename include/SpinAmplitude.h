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

#include "DataAccessor.h"
#include "DecayChannel.h"
#include "QuantumNumbers.h"
#include "ReportsInitialStateParticle.h"
#include "StaticDataAccessor.h"

#include <array>
#include <cstdlib>
#include <memory>

namespace yap {

class InitialStateParticle;

/// \class SpinAmplitude
/// \brief Abstract base class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup SpinAmplitude Spin Amplitudes

class SpinAmplitude : public StaticDataAccessor
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

    /// @}

    /// \return precalculated complex amplitude
    std::complex<double> amplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return Amplitude_->value(d, symmetrizationIndex(pc)); }

    /// access cached spin amplitude
    std::shared_ptr<ComplexCachedDataValue>& amplitude()
    { return Amplitude_; }

    /// access cached spin amplitude (const)
    const std::shared_ptr<ComplexCachedDataValue>& amplitude() const
    { return Amplitude_; }

    /// Add symmetrization indices for ParticleCombination.
    /// Must be overloaded in derived classes to both conditionally add ParticleCombination's
    /// and to alter and multiply them (to accomodate helicity, for example)
    virtual ParticleCombinationVector addSymmetrizationIndices(std::shared_ptr<ParticleCombination> pc) = 0;

    /// grant friend access to SpinAmplitudeCache to create SpinAmplitude's and set InitialStateParticle
    template <class spin_amplitude> friend class SpinAmplitudeCache;

protected:

    /// check equality
    virtual bool equals(const SpinAmplitude& other) const;

    /// Constructor
    /// declared private to ensure SpinAmplitude's are only created by a SpinAmplitudeCache
    /// \param intial quantum numbers of Initial-state
    /// \param final1 quantum numbers of first daughter
    /// \param final2 quantum numbers of second daughter
    /// \param l orbital angular momentum
    /// \param isp InitialStateParticle to which this SpinAmplitude belongs
    SpinAmplitude(const QuantumNumbers& initial,
                  const QuantumNumbers& final1,
                  const QuantumNumbers& final2,
                  unsigned l,
                  InitialStateParticle* isp);

private:

    /// Initial-state quantum numbers
    QuantumNumbers InitialQuantumNumbers_;

    /// array of final-state quantum numbers
    std::array<QuantumNumbers, 2> FinalQuantumNumbers_;

    /// orbital angular momentum
    unsigned L_;

    /// Cached complex spin amplitude
    std::shared_ptr<ComplexCachedDataValue> Amplitude_;

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

}

#endif
