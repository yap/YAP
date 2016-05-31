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

#include "fwd/Model.h"
#include "fwd/Particle.h"
#include "fwd/ParticleCombination.h"

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

    /// \return whether ParticleCombination is for a final state particle
    bool isFinalStateParticle() const
    { return Daughters_.empty() and Indices_.size() == 1; }

    /// \return top of decay tree this ParticleCombination belongs to
    std::shared_ptr<ParticleCombination> origin();

    /// \return vector of all leaves of decay tree below this ParticleCombination
    ParticleCombinationVector leaves();

    /// \return whether all leaves are final state particles
    bool decaysToFinalStateParticles() const;

    /// Checks consistency of object
    bool consistent() const;

    /// Create combinations.
    /// Example:
    /// P = { {f_0: [(01) -> (0) + (1)]; f_0: [(03) -> (0) + (3)]; f_0: [(21) -> (2) + (1)]; f_0: [(23) -> (2) + (3)]}
    ///       {pi+: [(0)]; pi+: [(2)]}
    ///       {pi-: [(1)]; pi-: [(3)]} }
    /// Result:
    /// { {nullptr: [(01) -> (0) + (1); (2); (3)]},
    ///   {nullptr: [(03) -> (0) + (3); (2); (1)]},
    ///   {nullptr: [(21) -> (2) + (1); (0); (3)]},
    ///   {nullptr: [(23) -> (2) + (3); (0); (1)]} }
    friend std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > combinations(std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > >& P, Model* m);

    /// Create combinations.
    /// Example:
    /// A = {pi+: [(0)]; pi+: [(2)]}
    /// B = { {pi-: [(1)]; pi-: [(3)]} }
    /// Result:
    /// { {nullptr: [(0); (1)]}
    ///   {nullptr: [(0); (3)]}
    ///   {nullptr: [(2); (1)]}
    ///   {nullptr: [(2); (3)]} }
    /// Results with overlapping indices will be omitted
    friend std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > combinations(std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> >& A, std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > B, Model* m);

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

    /// \name private constructors
    /// for valid use of shared_from_this()
    /// @{

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

    /// @}

    /// \name Equivalence-checking structs
    /// @{

public:

    /// \struct Equiv
    /// \brief base class for equivalence (with functor), compares shared_ptr's only
    struct Equiv {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
        { return A == B; }
    };

    /// \struct EquivByOrderedContent
    /// \brief Checks objects referenced by shared pointers, check indices only
    struct EquivByOrderedContent : Equiv {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivDown
    /// \brief Checks objects referenced by shared pointers,
    /// check self and all daughters (down the decay tree) for equality
    struct EquivDown : EquivByOrderedContent {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivUp
    /// \brief Check objects referenced by shared pointers,
    /// check self and parent (up the decay tree) for equality
    struct EquivUp : EquivByOrderedContent {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivUpAndDown
    /// \brief Check objects referenced by shared pointers,
    /// check self, all daughters (down-), and parent (up the decay tree) for equality
    struct EquivUpAndDown : EquivDown {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivByOrderlessContent
    /// \brief Check objects referenced by shared pointers,
    /// check indices only, disregarding order
    struct EquivByOrderlessContent : Equiv {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivDownByOrderlessContent
    /// \brief Check objects referenced by shared pointers,
    /// check indices only, disregarding order, and check daughters (but not daughter's daughters)
    /// Use e.g. for breakup momenta
    struct EquivDownByOrderlessContent : EquivByOrderlessContent {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivByReferenceFrame
    /// \brief Check objects referenced by shared pointers,
    /// Checks parents (and up) for orderless content
    /// Returns equivalent for content sitting in same reference frame.
    struct EquivByReferenceFrame : EquivByOrderlessContent {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };

    /// \struct EquivZemach

    /// \brief Check objects reference by shared pointers; treats all
    /// two-particle states as equal, compares three particle states
    /// by unordered content of resonance, and throws on 4 or more
    /// particles
    struct EquivZemach : Equiv {
        virtual bool operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const override;
    };


    /// \name Static Comparison objects
    static Equiv equivBySharedPointer;
    static EquivDown equivDown;
    static EquivUp equivUp;
    static EquivUpAndDown equivUpAndDown;
    static EquivByOrderedContent equivByOrderedContent;
    static EquivByOrderlessContent equivByOrderlessContent;
    static EquivDownByOrderlessContent equivDownByOrderlessContent;
    static EquivByReferenceFrame equivByReferenceFrame;
    static EquivZemach equivZemach;

/// @}

};

/// Get indices listed as string
std::string indices_string(const ParticleCombination& pc);

/// convert ParticleCombination to string
std::string to_string(const ParticleCombination& pc);

/// convert ParticleCombinationVector to string
std::string to_string(const ParticleCombinationVector& PCs);

/// convert vector<ParticleCombinationVector> to string
std::string to_string(const std::vector<ParticleCombinationVector>& PCs);

/// convert vector<vector<ParticleCombinationVector>> to string
std::string to_string(const std::vector<std::vector<ParticleCombinationVector>>& PCs);

/// convert pair<shared_ptr<Particle>, ParticleCombinationVector> to string
std::string to_string(const std::pair<std::shared_ptr<Particle>, ParticleCombinationVector>& pc);

/// convert vector<pair<shared_ptr<Particle>, ParticleCombinationVector> > to string
std::string to_string(const std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> >& PCs);

/// convert vector<vector<pair<shared_ptr<Particle>, ParticleCombinationVector> > > to string
std::string to_string(const std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > >& PCs);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombination& PC)
{ os << to_string(PC); return os; }

}

#endif
