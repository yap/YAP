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

#include <memory>
#include <set>
#include <string>

namespace yap {

/// \class ParticleCombinationCache
/// \brief Caches list of ParticleCombination's
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombinationCache
{
public:

    /// cache storage type
    using cache_type = std::set<std::weak_ptr<const ParticleCombination>, std::owner_less<std::weak_ptr<const ParticleCombination> > >;

    /// Default constructor
    ParticleCombinationCache() = default;

    /// Construct cache from vector of ISP ParticleCombination's
    ParticleCombinationCache(std::vector<std::shared_ptr<ParticleCombination> > ispPCs);

    /// check if cache contains element equating to pc
    /// \param pc shared_ptr to ParticleCombination to search for equivalent to
    cache_type::key_type find(std::shared_ptr<const ParticleCombination> pc) const;

    /// check if cache contains element equating to pc
    /// \param pc Pointer to ParticleCombination to search for equivalent to
    cache_type::key_type find(const ParticleCombination* pc) const
    { return find(std::shared_ptr<const ParticleCombination>(pc)); }

    /// check if cache contains element equating to pc
    /// \param pc ParticleCombination to search for equivalent to
    cache_type::key_type find(const ParticleCombination& pc) const;

    /// check if cache contains element equating to pc
    /// \param pc vector of ParticleIndex's to build ParticleCombination frrom for checking equivalence
    cache_type::key_type find(const std::vector<ParticleIndex>& I) const;

    /// \return shared_ptr to ParticleCombination from Cache, if it exists, otherwise adds it to cache.
    /// \param i ParticleIndex for FSP
    std::shared_ptr<const ParticleCombination> operator[](std::shared_ptr<const ParticleCombination> pc);

    /// \return shared_ptr to FSP ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache
    /// Constructs ParticleCombination from index.
    /// \param i ParticleIndex for FSP
    std::shared_ptr<const ParticleCombination> operator[](ParticleIndex i)
    { return operator[](std::make_shared<ParticleCombination>(i)); }

    /// \return shared_ptr to composite ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache.
    /// constructs PC from decay to final state particles in vector.
    /// \param I vector of ParticleIndex for FSP
    std::shared_ptr<const ParticleCombination> operator[](const std::vector<ParticleIndex>& I);

    /// \return shared_ptr to composite ParticleCombination from Cache, if it exists, otherwise constructs and adds it to cache.
    /// Constructs from ParticleCombinations for daughters
    /// \param c vector of shared_ptr's to ParticleCombination objects describing new ParticleCombination
    std::shared_ptr<const ParticleCombination> operator[](ParticleCombinationVector c)
    { return operator[](std::make_shared<ParticleCombination>(c)); }

    /// remove expired Cache_ elements
    void removeExpired();

    /// Check consistency of cache.
    /// Calls removeExpired()
    bool consistent() const;

    /// \name access to cache
    /// @{

    /// \return iterator to begin
    cache_type::iterator begin()
    { return Cache_.begin(); }

    /// \return const_iterator to begin
    cache_type::const_iterator begin() const
    { return Cache_.begin(); }

    /// \return iterator to end
    cache_type::iterator end()
    { return Cache_.end(); }

    /// \return const_iterator to end
    cache_type::const_iterator end() const
    { return Cache_.end(); }

    /// @}

protected:

    /// set lineage: copy each daughter, add pc as parent to copy,
    /// swap copy for daughter, and call setLineage on each daughter.
    void setLineage(std::shared_ptr<ParticleCombination> pc);

    /// set of weak pointers to ParticleCombination's
    cache_type Cache_;
};

/// convert to string
std::string to_string(const ParticleCombinationCache& C);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombinationCache& C)
{ os << to_string(C); return os; }

}

#endif
