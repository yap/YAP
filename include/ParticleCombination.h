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
#include <set>
#include <vector>

namespace yap {

class ParticleCombination;

/// \typedef ParticleCombinationSet
using ParticleCombinationSet = std::set<std::shared_ptr<const ParticleCombination> >;

/// \typedef ParticleCombinationVector
using ParticleCombinationVector = std::vector<std::shared_ptr<const ParticleCombination> >;

/// \typedef ParticleCombinationMap
template<typename T>
using ParticleCombinationMap = std::map<std::shared_ptr<const ParticleCombination>, T,
      std::owner_less<std::shared_ptr<const ParticleCombination> > >;

/// \class ParticleCombination
/// \brief Stores combinations of ParticleIndex types
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombination
{
public:

    /// Default constructor
    ParticleCombination();

    /// Final-state-particle constructor
    ParticleCombination(ParticleIndex index, char twoLambda = 0);

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

    /// Get vector of daughters (const)
    //ParticleCombinationVector daughters() const;

    /// get parent
    const ParticleCombination* parent() const
    { return Parent_; }

    /// get parent share_ptr
    const std::shared_ptr<const ParticleCombination> sharedParent() const;

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
    /// \return Success of action
    bool addDaughter(std::shared_ptr<const ParticleCombination> daughter);

    /// Checks consistency of combination
    /// by checking for absence of duplicate entries
    bool consistent() const;

    /// cast into string
    operator std::string() const;

    /// check if this and B share one or more ParticleIndex's
    bool sharesIndices(std::shared_ptr<const ParticleCombination> B) const;

    /// check if B is a subset of this
    /// Checks indices, but not parents or daughters
    bool isSubset(std::shared_ptr<const ParticleCombination> B) const;

    /// create new daughters with this PC as parent
    void setParents();

    /// add a particle combination as parent
    /// do not use. This function is used by makeParticleCombinationSetWithParents().
    void setParent(ParticleCombination* parent);

    /// set 2 * helicity
    void setTwoLambda(char twoLambda)
    { TwoLambda_ = twoLambda; }


    /// \name ParticleCombination friends
    /// @{

    /// equality operator
    friend bool operator==(const ParticleCombination& A, const ParticleCombination& B);

    /// inequality operator
    friend bool operator!=(const ParticleCombination& A, const ParticleCombination& B)
    { return !(A == B); }

    /// @}

protected:

    /// Parent of the particle combination.
    ParticleCombination* Parent_;
    ParticleCombinationVector Daughters_;
    std::vector<ParticleIndex> Indices_;
    /// 2 * Helicity
    char TwoLambda_;


/// \name Static methods for creating/retrieving ParticleCombination's
/// @{

// Following code is for managing unique shared pointers for particle
// combinations across all of YAP

public:

    /// return existing shared_ptr for final-state-particle ParticleCombination, if exists; otherwise creates and returns
    /// \param i ParticleIndex for FSP
    static std::shared_ptr<const ParticleCombination> uniqueSharedPtr(std::shared_ptr<const ParticleCombination> pc);

    /// return existing shared_ptr for final-state-particle ParticleCombination, if exists; otherwise creates and returns
    /// \param i ParticleIndex for FSP
    static std::shared_ptr<const ParticleCombination> uniqueSharedPtr(ParticleIndex i);

    /// return existing shared_ptr for final-state-particle ParticleCombination, if exists; otherwise creates and returns
    /// \param i ParticleIndex for FSP
    static std::shared_ptr<const ParticleCombination> uniqueSharedPtr(std::vector<ParticleIndex> I);

    /// return existing shared_ptr for ParticleCombination, if exists; otherwise creates and returns
    /// \param c vector of shared_ptr's to ParticleCombination objects describing new ParticleCombination
    static std::shared_ptr<const ParticleCombination> uniqueSharedPtr(ParticleCombinationVector c);

    /// return the particleCombination set
    static const ParticleCombinationSet& particleCombinationSet()
    { return ParticleCombinationSet_; }

    /// make a new particle combination set with parents set
    static void makeParticleCombinationSetWithParents(std::vector<std::shared_ptr<ParticleCombination> > initialStateParticleCombinations);

    static void printParticleCombinationSet();

private:

    /// \todo Move to ISP, make no longer static
    /// Static set of all particle combinations created throughout code
    static ParticleCombinationSet ParticleCombinationSet_;

/// @}

/// \name Comparison structs
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

/// convert ParticleCombination to string
std::string to_string(const ParticleCombination& pc);

}

#endif
