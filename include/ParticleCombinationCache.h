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

#ifndef yap_ParticleCombinationCache_h
#define yap_ParticleCombinationCache_h

#include "ParticleCombination.h"
#include "ParticleIndex.h"
#include "WeakPtrCache.h"

#include <memory>
#include <set>
#include <string>

namespace yap {

/// \class ParticleCombinationCache
/// \brief Caches list of ParticleCombination's
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombinationCache : public WeakPtrCache<const ParticleCombination>
{
public:

    /// implements equivalence checking
    bool equiv(const shared_ptr_type& A, const shared_ptr_type& B) const override
    { return ParticleCombination::equivUpAndDown(A, B); }

    /// Default constructor
    ParticleCombinationCache() = default;

    /// Construct cache from vector of ISP's ParticleCombination's
    ParticleCombinationCache(std::vector<shared_ptr_type> V);

    using WeakPtrCache::find;

    /// check if cache contains element equating to pc
    /// \param pc Pointer to ParticleCombination to search for equivalent to
    weak_ptr_type find(ParticleCombination* pc) const
    { return find(shared_ptr_type(pc)); }

    /// check if cache contains element equating to pc
    /// \param I vector of ParticleIndex's to build ParticleCombination from for checking equivalence
    weak_ptr_type find(const std::vector<ParticleIndex>& I) const;

    using WeakPtrCache::operator[];

    /// \return shared_ptr to FSP ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache
    /// Constructs ParticleCombination from index.
    /// \param i ParticleIndex for FSP
    shared_ptr_type operator[](ParticleIndex i)
    { return operator[](std::make_shared<type>(i)); }

    /// \return shared_ptr to composite ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache.
    /// constructs PC from decay to final state particles in vector.
    /// \param I vector of ParticleIndex for FSP
    shared_ptr_type operator[](const std::vector<ParticleIndex>& I);

    /// \return shared_ptr to composite ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache.
    /// Constructs from ParticleCombinations for daughters
    /// \param c vector of shared_ptr's to ParticleCombination objects describing new ParticleCombination
    shared_ptr_type operator[](ParticleCombinationVector c)
    { return operator[](std::make_shared<type>(c)); }

    /// Check consistency of cache.
    bool consistent() const;

protected:

    /// set lineage: copy each daughter, add pc as parent to copy,
    /// swap copy for daughter, and call setLineage on each daughter.
    void setLineage(shared_ptr_type pc);

    /// add to cache
    void addToCache(shared_ptr_type pc) override;

};

/// convert to string
std::string to_string(const ParticleCombinationCache& C);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombinationCache& C)
{ os << to_string(C); return os; }

}

#endif
