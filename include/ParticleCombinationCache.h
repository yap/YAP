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

class ParticleCombinationCache : public WeakPtrCache<ParticleCombination>
{
public:

    /// implements equivalence checking
    bool equiv(const shared_ptr_type& A, const shared_ptr_type& B) const override
    { return ParticleCombination::equivUpAndDown(A, B); }

    /// Default constructor
    ParticleCombinationCache() = default;

    /// Construct cache from vector of ISP's ParticleCombination's
    ParticleCombinationCache(std::vector<shared_ptr_type> V);

    /// retrieve or create final-state particle ParticleCombination
    /// \param index Index of particle
    shared_ptr_type fsp(ParticleIndex index)
    { return operator[](create_fsp(index)); }

    /// retrieve or create composite particle from daughters.
    /// copies daughters into composite, setting copies' parents = shared_from_this()
    /// \param D ParticleCombinationVector of daughters to create composite from
    shared_ptr_type composite(const ParticleCombinationVector& D)
    { return operator[](create_composite(D)); }

    using WeakPtrCache::find;

    /*
        /// retrieve final-state particle ParticleCombination
        /// Does not add to the cache if ParticleCombination is not found.
        /// \param index Index of particle
        weak_ptr_type find(ParticleIndex index)
        { return find(create_fsp(index)); }

        /// retrieve copy of ParticleCombination with new spin projection
        /// Does not add to the cache if ParticleCombination is not found.
        /// \param other ParticleCombination to copy
        weak_ptr_type find(const ParticleCombination& other)
        { return find(create_copy(other)); }
    */

    /// retrieve composite particle ParticleCombination from cache.
    /// Does not add to the cache if ParticleCombination is not found.
    /// \param D vector of daughters to construct ParticleCombination from
    /// \return weak_ptr to ParticleCombination; is empty if not found.
    weak_ptr_type find(const ParticleCombinationVector& D)
    { return find(create_composite(D)); }

    /// retrieves first entry matching vector of particle indices by unordered content
    /// Does not add to the cache if a match is not found.
    /// \param I vector of ParticleIndex's to look for
    /// \return weak_ptr to ParticleCombination; is empty if not found.
    weak_ptr_type findByUnorderedContent(const std::vector<ParticleIndex>& I) const;

    /// Check consistency of cache.
    bool consistent() const;

    /// stream the cache elements as a table
    virtual std::ostream& print(std::ostream& os) const override;

protected:

    /* /// set lineage: copy each daughter, add pc as parent to copy, */
    /* /// swap copy for daughter, and call setLineage on each daughter. */
    /* void setLineage(shared_ptr_type pc); */

private:

    /// add to cache
    void addToCache(shared_ptr_type pc) override;

    /// create final-state particle, it is not yet added to cache
    shared_ptr_type create_fsp(ParticleIndex index) const
    { return shared_ptr_type(new ParticleCombination(index)); }

    /// create copy, it is not yet added to cache
    shared_ptr_type create_copy(const ParticleCombination& other) const;

    /// create composite ParticleCombination, it is not yet added to cache
    shared_ptr_type create_composite(const ParticleCombinationVector& D) const;

};

}

#endif
