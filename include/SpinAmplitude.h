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

#include "fwd/SpinAmplitude.h"

#include "fwd/ParticleCombination.h"
#include "fwd/Spin.h"
#include "fwd/StatusManager.h"

#include "CachedValue.h"
#include "MathUtilities.h"
#include "StaticDataAccessor.h"

#include <array>
#include <complex>
#include <iostream>
#include <memory>

namespace yap {

/// \class SpinAmplitude
/// \brief Abstract base class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup SpinAmplitude Spin Amplitudes
class SpinAmplitude : public StaticDataAccessor
{
public:

    /// \typedef AmplitudeSubmap
    /// \brief maps SpinProjectionVector to ComplexCachedValue
    using AmplitudeSubmap = std::map<SpinProjectionVector, std::shared_ptr<ComplexCachedValue> >;

    /// \typedef AmplitudeMap
    /// \brief maps parent spin projectin to AmplitudeSubmap
    using AmplitudeMap = std::map<int, AmplitudeSubmap>;

    /// cast into string
    virtual operator std::string() const;

    /// \name Getters
    /// @{

    /// \return initial spin * 2 (const)
    const unsigned initialTwoJ() const
    { return InitialTwoJ_; }

    /// \return SpinVector of daughters (const)
    const SpinVector& finalTwoJ() const
    { return FinalTwoJ_; }

    /// \return orbital angular momentum
    const unsigned L() const
    { return L_; }

    /// \return total spin angular momentum */
    const unsigned twoS() const
    { return TwoS_; }

    /// \return set of (twice the) spin projections of initial state
    const std::set<int> twoM() const;

    /// \return set of daughter spin projections for given initial state spin projection
    const std::set<SpinProjectionVector> twoM(int two_M) const;

    /// @}

    /// Calculate spin amplitude for caching.
    /// Must be overrided in derived classes.
    /// \param two_M 2 * spin projection of parent
    /// \param two_m SpinProjectionVector of daughters
    /// \param d DataPoint to retrieve data from for calculation
    /// \param pc ParticleCombination to calculate for
    virtual const std::complex<double> calc(int two_M, const SpinProjectionVector& two_m,
                                            const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const = 0;

    /// Loops over particle combinations (pc) and all (M, vec{m}) combinations
    /// and call calc(M, \vec{m}, d, pc) when necessary
    /// \param d DataPoint to calculate into
    /// \param sm StatusManager to update
    void calculate(DataPoint& d, StatusManager& sm) const override;

    /// \return precalculated complex amplitude
    /// \param d DataPoint to retrieve value from
    /// \param pc ParticleCombination to retrieve value for
    /// \param two_M 2 * spin projection of parent
    /// \param two_m SpinProjectionVector of daughters
    virtual const std::complex<double> amplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                         int two_M, const SpinProjectionVector& two_m) const
    { return Amplitudes_.at(two_M).at(two_m)->value(d, symmetrizationIndex(pc)); }

    /// \return a string naming the formalism used for the SpinAmplitude calculation
    virtual std::string formalism() const = 0;

    /// grant friend access to SpinAmplitudeCache to set Model
    friend class SpinAmplitudeCache;

    /// grant friend access to DecayChannel to call addParticleCombination
    friend class DecayChannel;

protected:

    /// add spin amplitude for transition from state with initial
    /// projection to states with final projections.
    /// \param two_M twice the spin projection of the initial state
    /// \param two_m SpinProjectionVector of daughters
    /// \param addNULLCachedValue set true for use from UnitSpinAmplitude
    virtual void addAmplitude(int two_M, const SpinProjectionVector& two_m, bool addNULLCachedValue = false);

    /// check equivalence: only check spins and angular momenta
    virtual bool equalTo(const SpinAmplitude& other) const;

    /// check equality: calls #equalTo and checks symmetrizationIndices
    virtual bool equals(const SpinAmplitude& other) const;

    /// Constructor
    /// declared private to ensure SpinAmplitude's are only created by a SpinAmplitudeCache
    /// \param two_J twice the spin of Initial-state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    /// \param equal ParticleCombination equality struct for determining index assignments
    SpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s,
                  const ParticleCombinationEqualTo& equal);

private:

    /// Initial-state spin * 2
    unsigned InitialTwoJ_;

    /// SpinVector of daughters
    SpinVector FinalTwoJ_;

    /// orbital angular momentum
    unsigned L_;

    /// twice the total spin angular momentum
    unsigned TwoS_;

    /// Cached complex spin amplitude
    AmplitudeMap Amplitudes_;

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

/// convert to string
std::string to_string(const SpinAmplitudeVector& saV);

}

#endif
