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

#include "ParticleIndex.h"

#include <map>
#include <memory>
#include <ostream>
#include <vector>

namespace yap {

class ParticleCombination;

/// \typedef ParticleCombinationVector
using ParticleCombinationVector = std::vector<std::shared_ptr<const ParticleCombination> >;

/// \typedef ParticleCombinationMap
template<typename T>
using ParticleCombinationMap = std::map<std::shared_ptr<const ParticleCombination>, T,
      std::owner_less<std::shared_ptr<const ParticleCombination> > >;

/// \class ParticleCombination
/// \brief Stores combinations of ParticleIndex types
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombination // : public std::enable_shared_from_this<ParticleCombination>
{
public:

    /// Default constructor
    ParticleCombination() = default;

    /// Final-state-particle constructor
    ParticleCombination(ParticleIndex index, char twoLambda = 0) : Indices_(1, index), TwoLambda_(twoLambda) {}

    /// Resonance particle constructor
    ParticleCombination(ParticleCombinationVector c, char twoLambda = 0);

    /// \name Getters
    /// @{

    /// Get vector of indices
    const std::vector<ParticleIndex>& indices() const
    { return Indices_; }

    /// Get vector of daughters
    const ParticleCombinationVector& daughters() const
    { return Daughters_; }

    /// get parent
    std::shared_ptr<const ParticleCombination> parent() const
    { return std::const_pointer_cast<const ParticleCombination>(Parent_.lock()); }

    /// get 2 * helicity
    const char twoLambda() const
    { return TwoLambda_; }

    /// @}

    /// \name Get info on type
    /// @{

    bool isFinalStateParticle() const
    { return Daughters_.empty() and Indices_.size() == 1; }

    /// @}

    /// Add daughter ParticleCombination
    /// \param daughter Shared pointer to ParticleCombination object representing a daughter
    void addDaughter(std::shared_ptr<const ParticleCombination> daughter);

    /// Checks consistency of combination
    /// by checking for absence of duplicate entries
    bool consistent() const;

    /// set 2 * helicity
    void setTwoLambda(char twoLambda)
    { TwoLambda_ = twoLambda; }

    /// grant friend access for setting lineage
    friend class ParticleCombinationCache;

protected:

    /// Parent of the particle combination.
    std::weak_ptr<ParticleCombination> Parent_;

    /// vector of daughters
    ParticleCombinationVector Daughters_;

    /// vector indices of daughters
    std::vector<ParticleIndex> Indices_;

    /// 2 * Helicity
    char TwoLambda_;

    /// raw pointer to ParticleCombinationCache managing this object
    /// \todo Do we need this? Do we really need to check if the PC is in the cache during consistency checking?
    ParticleCombinationCache Cache_(nullptr);

    /// \name Equivalence-checking structs
    /// @{

public:

    /// \struct Equiv
    /// \brief base class for equivalence (with functor), compares shared_ptr's only
    struct Equiv {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
        { return A == B; }
    };

    /// \struct EquivByOrderedContent
    /// \brief Checks objects referenced by shared pointers, check indices only
    /// Does NOT compare helicity
    struct EquivByOrderedContent : Equiv {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivDownButLambda
    /// \brief Checks objects referenced by shared pointers,
    /// check self and all daughters (down the decay tree) for equality
    /// Does NOT compare helicity
    struct EquivDownButLambda : EquivByOrderedContent {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivDown
    /// \brief Checks objects referenced by shared pointers,
    /// check self and all daughters (down the decay tree) for equality
    /// Also compares helicity
    struct EquivDown : EquivByOrderedContent {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivUpButLambda
    /// \brief Check objects referenced by shared pointers,
    /// check self, and parents (up the decay tree) for equality
    /// Does NOT compare helicity
    struct EquivUpButLambda : EquivByOrderedContent {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivUpAndDownButLambda
    /// \brief Check objects referenced by shared pointers,
    /// check self, all daughters (down-), and parents (up the decay tree) for equality
    /// Does NOT compare helicity
    struct EquivUpAndDownButLambda : EquivDownButLambda {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivUpAndDown
    /// \brief Check objects referenced by shared pointers,
    /// check self, all daughters (down-), and parent (up the decay tree) for equality
    /// Also compares helicity
    struct EquivUpAndDown : EquivDown {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivByOrderlessContent
    /// \brief Check objects referenced by shared pointers,
    /// check indices only, disregarding order
    /// Does NOT compare helicity
    struct EquivByOrderlessContent : Equiv {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivDownByOrderlessContent
    /// \brief Check objects referenced by shared pointers,
    /// check indices only, disregarding order, and check daughters (but not daughter's daughters)
    /// Does NOT compare helicity
    /// Use e.g. for breakup momenta
    struct EquivDownByOrderlessContent : EquivByOrderlessContent {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \struct EquivByReferenceFrame
    /// \brief Check objects referenced by shared pointers,
    /// Checks parents (and up) for orderless content
    /// Does NOT compare helicity. Returns equivalent for content sitting in same reference frame.
    struct EquivByReferenceFrame : EquivByOrderlessContent {
        virtual bool operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const override;
    };

    /// \name Static Comparison objects
    static Equiv equivBySharedPointer;
    static EquivDown equivDown;
    static EquivDownButLambda equivDownButLambda;
    static EquivUpAndDown equivUpAndDown;
    static EquivUpButLambda equivUpButLambda;
    static EquivUpAndDownButLambda equivUpAndDownButLambda;
    static EquivByOrderedContent equivByOrderedContent;
    static EquivByOrderlessContent equivByOrderlessContent;
    static EquivDownByOrderlessContent equivDownByOrderlessContent;
    static EquivByReferenceFrame equivByReferenceFrame;

/// @}

};

/// equality operator
inline bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{ return ParticleCombination::equivUpAndDown(std::make_shared<ParticleCombination>(A), std::make_shared<ParticleCombination>(B)); }

/// inequality operator
inline bool operator!=(const ParticleCombination& A, const ParticleCombination& B)
{ return !(A == B); }

/// convert ParticleCombination to string
std::string to_string(const ParticleCombination& pc);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombination& PC)
{ os << to_string(PC); return os; }

}

#endif
