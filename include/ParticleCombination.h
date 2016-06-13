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

#include <algorithm>
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

/// \return if the given ParticleCombination is in ParticleCombinations
/// \param PCs ParticleCombinations
/// \param c ParticleCombination to look for equivalent of in ParticleCombinations
/// \param equiv ParticleCombination::Equiv object for checking equivalence
inline bool hasParticleCombination(const ParticleCombinationVector PCs,
                                   const std::shared_ptr<ParticleCombination>& c,
                                   const ParticleCombination::Equiv& equiv = ParticleCombination::equivBySharedPointer)
{ return (std::find_if(PCs.begin(), PCs.end(), [&](const ParticleCombinationVector::value_type & pc) { return equiv(pc, c); } )) != PCs.end(); }

/// Get indices listed as string
std::string indices_string(const ParticleCombination& pc);

/// convert ParticleCombination to string
std::string to_string(const ParticleCombination& pc);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombination& PC)
{ os << to_string(PC); return os; }

}

#endif
