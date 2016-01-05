#include "ParticleCombination.h"

#include "container_utils.h"
#include "Exceptions.h"
#include "logging.h"
#include "QuantumNumbers.h"

#include <algorithm>
#include <set>

namespace yap {

//-------------------------
void ParticleCombination::addDaughter(std::shared_ptr<const ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        FLOG(ERROR) << "daughter contains no indices.";
        throw exceptions::Exception("daughter contains no indices", "ParticleCombination::addDaughter");
    }

    /// Check that new daughter does not share content with other daughters?
    if (overlap(daughter->indices(), Indices_)) {
        FLOG(ERROR) << "daughter contains indices that are already in parent.";
        throw exceptions::Exception("daughter overlaps with other daughters", "ParticleCombination::addDaughter");
    }

    /// Check that new daughter doesn't already have parent set
    if (daughter->parent()) {
        FLOG(ERROR) << "daughter's parent is already set.";
        throw exceptions::Exception("daughter's parent is already set", "ParticleCombination::addDaughter");
    }

    // copy daughter
    auto d = new ParticleCombination(*daughter);
    
    // set daughter's parent to shared_from_this
    d->Parent_ = shared_from_this();
    
    // add daughter to vector
    Daughters_.emplace_back(d);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), Daughters_.back()->indices().begin(), Daughters_.back()->indices().end());
}

//-------------------------
bool ParticleCombination::consistent() const
{
    bool C = true;

    // should have no daughters or 2 or more daughters:
    if (Daughters_.size() == 1) {
        FLOG(ERROR) << "has only one daughter.";
        C &= false;
    }

    if (Indices_.empty()) {
        FLOG(ERROR) << "has no indices.";
        C &= false;
    }

    // if has daughters
    if (!Daughters_.empty()) {
        // Check Indices_ doesn't have duplicates
        // create unique_copy of Indices_ (as set)
        std::set<ParticleIndex> U(Indices_.begin(), Indices_.end());
        // check unique_copy-object's size == this object's size
        if (U.size() != Indices_.size()) {
            FLOG(ERROR) << "index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
            C &= false;
        }

        // count number of daughters with parent not set to this
        auto n = std::count_if(Daughters_.begin(), Daughters_.end(),
        [&](const ParticleCombinationVector::value_type & d) {return d->parent().get() != this;});
        if (n != 0) {
            FLOG(ERROR) << n << " daughters' parent not set to this ParticleCombination.";
            C &= false;
        }

        // check consistency of daughters
        std::for_each(Daughters_.begin(), Daughters_.end(),
        [&](const ParticleCombinationVector::value_type & d) {C &= d->consistent();});
    }
    // if Daugthers_ empty, should have one and only index (as FSP)
    else if (Indices_.size() != 1) {
        FLOG(ERROR) << "contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
        C &= false;
    }

    return C;
}

//-------------------------
std::string to_string(const ParticleCombination& pc)
{
    if (pc.indices().empty())
        return "(empty)";
    std::string s = "(";
    for (auto i : pc.indices())
        s += std::to_string(i);
    s += ", λ = " + spin_to_string(pc.twoLambda()) + ")";

    if (pc.daughters().empty() or (pc.daughters().size() == 2 and pc.indices().size() == 2))
        return s;

    s += " -> ";

    for (auto& d : pc.daughters()) {
        if (d->daughters().empty() or (d->daughters().size() == 2 and d->indices().size() == 2))
            s += to_string(*d) + ", ";
        else
            s += "[" + to_string(*d) + "], ";
    }
    s.erase(s.size() - 2, 2);
    return s;
}

//-------------------------
// Comparison stuff:

ParticleCombination::Equiv ParticleCombination::equivBySharedPointer;
ParticleCombination::EquivDown ParticleCombination::equivDown;
ParticleCombination::EquivDownButLambda ParticleCombination::equivDownButLambda;
ParticleCombination::EquivUpAndDown ParticleCombination::equivUpAndDown;
ParticleCombination::EquivUpButLambda ParticleCombination::equivUpButLambda;
ParticleCombination::EquivUpAndDownButLambda ParticleCombination::equivUpAndDownButLambda;
ParticleCombination::EquivByOrderedContent ParticleCombination::equivByOrderedContent;
ParticleCombination::EquivByOrderlessContent ParticleCombination::equivByOrderlessContent;
ParticleCombination::EquivDownByOrderlessContent ParticleCombination::equivDownByOrderlessContent;

//-------------------------
bool ParticleCombination::EquivByOrderedContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // Check indices
    if (A->indices().size() != B->indices().size())
        return false;
    if (!std::equal(A->indices().begin(), A->indices().end(), B->indices().begin()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivDownButLambda::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        if (!operator()(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivDown::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    /// check lambda
    if (A->TwoLambda_ != B->TwoLambda_)
        return false;

    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        if (!operator()(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpButLambda::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // check parent
    if (! ParticleCombination::equivUpButLambda(A->parent(), B->parent()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpAndDownButLambda::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivDownButLambda(A, B))
        return false;

    // check parent
    if (! ParticleCombination::equivUpButLambda(A->parent(), B->parent()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpAndDown::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivDown(A, B))
        return false;

    // check parent
    if (A->parent() != B->parent())
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivByOrderlessContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check size of indices vectors
    if (A->indices().size() != B->indices().size())
        return false;

    // check contents of indices vectors
    // (creating a set will sort entries for easy comparison,
    // since order doesn't matter)
    std::set<ParticleIndex> a(A->indices().begin(), A->indices().end());
    std::set<ParticleIndex> b(B->indices().begin(), B->indices().end());

    return std::equal(a.begin(), a.end(), b.begin());
}

//-------------------------
bool ParticleCombination::EquivDownByOrderlessContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderlessContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }

    unsigned matches(0);
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        for (unsigned j = 0; j < B->daughters().size(); ++j)
            if (ParticleCombination::equivByOrderlessContent(A->daughters()[i], B->daughters()[j])) {
                ++matches;
                break;
            }

    if (matches == A->daughters().size())
        // a match!
        return true;

    return false;
}

//-------------------------
bool ParticleCombination::EquivByReferenceFrame::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    // if both are nullptr, also return true
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderlessContent(A->parent(), B->parent()))
        return false;

    return operator()(A->parent(), B->parent());
}

}
