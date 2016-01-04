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

#include "SpinAmplitude.h"
#include "WeakPtrCache.h"

#include <memory>

namespace yap {

/// \class SpinAmplitudeCache
/// \brief Caches SpinAmplitudes
/// \author Johannes Rauch, Daniel Greenwald

class SpinAmplitudeCache : public WeakPtrCache<SpinAmplitude>
{
public:

    /// equivalence
    bool equiv(const std::shared_ptr<SpinAmplitude>& A, const std::shared_ptr<SpinAmplitude>& B) const override
    { return (A.get() == B.get()) or (*A == *B); }

    /// Check consistency of cache. Skips expired entries.
    bool consistent() const
    {
        bool C = true;
        for (auto it = begin(); it != end(); ++it)
            if (!it->expired())
                C &= it->lock()->consistent();
        return C;
    }

};

/*
/// convert to string
std::string to_string(const SpinAmplitudeCache& C);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const SpinAmplitudeCache& C)
{ os << to_string(C); return os; }
*/

}

#endif
