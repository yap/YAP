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

#ifndef yap_SpinAmplitudeCache_h
#define yap_SpinAmplitudeCache_h

#include "ReportsInitialStateParticle.h"
#include "SpinAmplitude.h"
#include "WeakPtrCache.h"

#include <memory>

namespace yap {

class InitialStateParticle;

/// \class SpinAmplitudeCache
/// \brief Caches SpinAmplitudes
/// \author Johannes Rauch, Daniel Greenwald
///
/// Templating here insures that all SpinAmplitude's created for an
/// InitialStateParticle have the same formalism
///
/// \tparam spin_amplitude Class for constructing SpinAmplitude's from
template <class spin_amplitude>
class SpinAmplitudeCache :
    public WeakPtrCache<SpinAmplitude>,
    public ReportsInitialStateParticle
{
public:

    /// Constructor
    /// \param isp raw pointer to InitialStateParticle this cache belongs to
    SpinAmplitudeCache(InitialStateParticle* isp = nullptr) :
        WeakPtrCache(), InitialStateParticle_(isp) {}

    /// equivalence
    bool equiv(const std::shared_ptr<SpinAmplitude>& A, const std::shared_ptr<SpinAmplitude>& B) const override
    { return (A.get() == B.get()) or A->equiv(*B); }

    /// retrieve or create SpinAmplitude
    /// \param intial quantum numbers of Initial-state
    /// \param final1 quantum numbers of first daughter
    /// \param final2 quantum numbers of second daughter
    /// \param l orbital angular momentum
    /// \param two_s 2 * the total spin angular momentum
    std::shared_ptr<SpinAmplitude> spinAmplitude(const QuantumNumbers& initial,
            const QuantumNumbers& final1,
            const QuantumNumbers& final2,
            unsigned l, unsigned two_s)
    {
        auto retVal = operator[](std::shared_ptr<spin_amplitude>(new spin_amplitude(initial, final1, final2, l, two_s, initialStateParticle())));
        retVal->setInitialStateParticle(initialStateParticle());
        return retVal;
    }

    /// Check consistency of cache. Skips expired entries.
    bool consistent() const
    {
        bool C = true;
        for (auto it = begin(); it != end(); ++it)
            if (!it->expired())
                C &= it->lock()->consistent();
        return C;
    }

    /// \return raw pointer to owning InitialStateParticle
    InitialStateParticle* initialStateParticle() override
    { return InitialStateParticle_; }

private:

    /// raw pointer to InitialStateParticle this cache belongs to
    InitialStateParticle* InitialStateParticle_;

};

}

#endif
