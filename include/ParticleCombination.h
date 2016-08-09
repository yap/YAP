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

#ifndef yap_ParticleCombination_h
#define yap_ParticleCombination_h

#include "fwd/ParticleCombination.h"

#include "fwd/Model.h"

#include <algorithm>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace yap {

/// \class ParticleCombination
/// \brief Stores combinations of particle indices
/// \author Johannes Rauch, Daniel Greenwald
///
/// Constructors are private. New ParticleCombination objects are
/// created through the ParticleCombinationCache
class ParticleCombination : public std::enable_shared_from_this<ParticleCombination>
{
private:

    /// default constructor
    ParticleCombination() = default;

    /// Final-state-particle constructor, see ParticleCombinationCache::fsp for details
    ParticleCombination(unsigned index) : Indices_(1, index) {}

    /// Copy constructor is deleted
    ParticleCombination(const ParticleCombination&) = delete;

    /// Move constructor is deleted
    ParticleCombination(ParticleCombination&&) = delete;

    /// Copy assignment is deleted
    ParticleCombination& operator=(const ParticleCombination&) = delete;

    /// Move assignment is deleted
    ParticleCombination& operator=(ParticleCombination&&) = delete;

public:

    /// \name Getters
    /// @{

    /// Get vector of indices (const)
    const std::vector<unsigned>& indices() const
    { return Indices_; }

    /// Get vector of daughters (const)
    const ParticleCombinationVector& daughters() const
    { return Daughters_; }

    /// get parent (const)
    std::shared_ptr<ParticleCombination> parent() const
    { return Parent_.lock(); }

    /// @}

    /// Checks consistency of object
    bool consistent() const;

    /// grant friend access to ParticleCombinationCache for creating ParticleCombination's
    friend class ParticleCombinationCache;

protected:

    /// Add daughter ParticleCombination
    /// \param daughter Shared pointer to ParticleCombination object representing a daughter
    void addDaughter(std::shared_ptr<ParticleCombination> daughter);

private:

    /// Parent of the particle combination.
    std::weak_ptr<ParticleCombination> Parent_;

    /// vector of daughters
    ParticleCombinationVector Daughters_;

    /// vector indices of daughters
    std::vector<unsigned> Indices_;

};

/// \name ParticleCombinationEqualTo functions
/// @{

/// compare shared_ptr's to ParticleCombination, returning true always
constexpr bool equal_always(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{ return true; }

/// compare shared_ptr's to ParticleCombination by shared_ptr only
inline bool equal_by_shared_pointer(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{ return A == B; }

/// compare shared_ptr's to ParticleCombination by indices only
bool equal_by_ordered_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// compare shared_ptr's to ParticleCombination by indices, disregarding order
bool equal_by_orderless_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// compare shared_ptr's to ParticleCombination by checking selves and then all daughters
bool equal_down(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// compare shared_ptr's to ParticleCombination by checking selves and then parents
bool equal_up(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// compare shared_ptr's to ParticleCombination by checking selves, daughters, and parents
bool equal_up_and_down(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// compare shared_ptr's to ParticleCombination by checking selves
/// with #equal_by_orderless_content and then checking (only) one generation down
bool equal_down_by_orderless_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B);

/// @}

/// \return whether all members of ParticleCombinationVector are non-overlapping with each other
/// \param pcv Vector check in
bool disjoint(const ParticleCombinationVector& pcv);

/// \return whether ParticleCombination is for a final state particle
inline const bool is_final_state_particle_combination(const ParticleCombination& pc)
{ return pc.daughters().empty() and pc.indices().size() == 1; }

/// \return top of decay tree this ParticleCombination belongs to
inline const ParticleCombination& origin(const ParticleCombination& pc)
{ return pc.parent() ? origin(*pc.parent()) : pc; }

/// \return whether pc is a pc of an initial state particle
const bool is_initial_state_particle_combination(const ParticleCombination& pc, const Model& m);

/// \return whether pc is or is from a pc of an initial state particle
inline const bool is_from_initial_state_particle_combination(const ParticleCombination& pc, const Model& m)
{ return is_initial_state_particle_combination(origin(pc), m); }

/// only keep particleCombinations with the highest number of indices in their top-most parent
void prune_particle_combinations(ParticleCombinationSet& PCs);

/// Get indices listed as string
std::string indices_string(const ParticleCombination& pc, std::string before = "(", std::string after = ")");

/// convert ParticleCombination to string
std::string to_string(const ParticleCombination& pc);

/// convert ParticleCombination with top-most parent to string
std::string to_string_with_parent(const ParticleCombination& pc);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombination& PC)
{ os << to_string(PC); return os; }

}

#endif
