#include "ParticleCombination.h"

#include "logging.h"
#include "MathUtilities.h"
#include "SpinUtilities.h"

#include <algorithm>
#include <assert.h>

namespace yap {

//-------------------------
ParticleCombination::ParticleCombination() :
    Parent_(nullptr)
{
}

//-------------------------
ParticleCombination::ParticleCombination(ParticleIndex index, char twoLambda) :
    Parent_(nullptr),
    Indices_(1, index),
    TwoLambda_(twoLambda)
{
}

//-------------------------
ParticleCombination::ParticleCombination(ParticleCombinationVector c, char twoLambda) :
    Parent_(nullptr),
    TwoLambda_(twoLambda)
{
    for (auto& d : c)
        addDaughter(d);
}

//-------------------------
std::string to_string(const ParticleCombination& pc)
{
    if (pc.indices().empty())
        return "(empty)";
    std::string s = "(";
    std::for_each(pc.indices().begin(), pc.indices().end(), [&](const ParticleIndex & i) {s += std::to_string(i) + ", ";});
    if (!pc.indices().empty())
        s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}

//-------------------------
void ParticleCombination::addDaughter(std::shared_ptr<const ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        FLOG(ERROR) << "daughter contains no indices.";
        throw std::invalid_argument("daughter");
    }

    /// Check that new daughter does not share content with other daughters?
    if (overlap(daughter->indices(), Indices_)) {
        FLOG(ERROR) << "daughter contains indices that are already in parent.";
        throw std::invalid_argument("daughter");
    }

    // add daughter to vector
    Daughters_.push_back(daughter);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), daughter->indices().begin(), daughter->indices().end());
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

    // if Daugthers_ empty, should have one and only index
    if (Daughters_.empty()) {
        if (Indices_.size() != 1) {
            FLOG(ERROR) << "contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
            C &= false;
        } else
            // don't need to check indices below
            return true;
    } else {

        // Check Indices_ doesn't have duplicates
        // create unique_copy of Indices_ (as set)
        std::set<ParticleIndex> U(Indices_.begin(), Indices_.end());
        // check unique_copy-object's size == this object's size
        if (U.size() != Indices_.size()) {
            FLOG(ERROR) << "index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
            C &= false;
        }

        // check if in ParticleCombinationCache_
        /// \todo make check_cached function
        if (std::find_if(ParticleCombinationCache_.begin(), ParticleCombinationCache_.end(),
                         [&](const ParticleCombinationCache::key_type& pc){return pc.lock().get() == this;}) == ParticleCombinationCache_.end()) {
            FLOG(ERROR) << "ParticleCombination is not in ParticleCombinationSet: " << std::string(*this)
                        << ((parent()) ? std::string(" from decay ") + std::string(*parent()) : " (no parent)");
            C & = false;
        }

        // count number of daughters with parent not set to this
        auto n = std::count_if(Daughters_.begin(), Daughters_.end(),
                               [&](const ParticleCombinationVector::value_type& d){return d->parent().get() != this;});
        if (n != 0) {
            FLOG(ERROR) << n << " daughters' parent not set to this ParticleCombination.";
            C &= false;
        }

        // check consistency of daughters
        std::for_each(Daughers_.begin(), Daughters.end(),
                      [&](const ParticleCombinationVector::value_type& d){C &= d->consistent();});
    }
    
    return C;
}

//-------------------------
ParticleCombination::operator std::string() const
{
    std::string result = "(";
    for (ParticleIndex i : Indices_)
        result += std::to_string(static_cast<unsigned>(i));
    /*std::ostringstream address;
    address << ", " << Parent_ << "->" << this;
    result += address.str();*/
    result += " λ=" + spinToString(TwoLambda_);
    result += ")";

    if (Daughters_.empty() || (Daughters_.size() == 2 and Indices_.size() == 2))
        return result;

    result += " -> ";

    for (auto& d : Daughters_) {
        if (d->daughters().empty() || (d->daughters().size() == 2 and d->indices().size() == 2))
            result += std::string(*d);
        else
            result += "[" + std::string(*d) + "]";

    }


    return result;
}

//-------------------------
bool ParticleCombination::sharesIndices(std::shared_ptr<const ParticleCombination> B) const
{
    for (ParticleIndex a : Indices_)
        for (ParticleIndex b : B->indices())
            if (a == b)
                return true;

    return false;
}

//-------------------------
bool ParticleCombination::isSubset(std::shared_ptr<const ParticleCombination> B) const
{
    for (ParticleIndex b : B->indices()) {
        bool found(false);
        for (ParticleIndex a : Indices_) {
            if (a == b) {
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }

    return true;
}

//-------------------------
void ParticleCombination::setParents()
{
    for (auto& daughter : Daughters_) {
        // make copy, set this as parent and get unique shared_ptr
        std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*daughter));
        copy->Parent_ = this;

        // call recursively
        copy->setParents();

        std::shared_ptr<const ParticleCombination> uniqueCopy = uniqueSharedPtr(copy);
        daughter.swap(uniqueCopy);
    }
}

//-------------------------
bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{
    return ParticleCombination::equivUpAndDown(std::make_shared<ParticleCombination>(A), std::make_shared<ParticleCombination>(B));
}

/////////////////////////
// Static stuff:

ParticleCombinationSet ParticleCombination::ParticleCombinationSet_;

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombination::uniqueSharedPtr(std::shared_ptr<const ParticleCombination> pc)
{
    for (auto& d : ParticleCombinationSet_)
        if (ParticleCombination::equivUpAndDown(pc, d))
            return d;

    ParticleCombinationSet_.insert(pc);
    return pc;
}

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombination::uniqueSharedPtr(ParticleIndex i)
{
    return uniqueSharedPtr(std::make_shared<ParticleCombination>(i));
}

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombination::uniqueSharedPtr(std::vector<ParticleIndex> I)
{
    ParticleCombinationVector V;
    for (ParticleIndex i : I)
        V.push_back(uniqueSharedPtr(i));
    return uniqueSharedPtr(V);
}

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombination::uniqueSharedPtr(ParticleCombinationVector c)
{
    return uniqueSharedPtr(std::make_shared<yap::ParticleCombination>(c));
}

//-------------------------
void ParticleCombination::makeParticleCombinationSetWithParents(std::vector<std::shared_ptr<ParticleCombination> > initialStateParticleCombinations)
{
    ParticleCombinationSet_.clear();

    for (auto& pc : initialStateParticleCombinations) {
        pc->setParents();
    }

    for (auto& pc : initialStateParticleCombinations)
        ParticleCombinationSet_.insert(pc);

    // check consistency
    for (auto& pc : ParticleCombinationSet_) {
        assert(pc->consistent());
    }
}

//-------------------------
void ParticleCombination::printParticleCombinationSet()
{
    std::cout << "ParticleCombination set:\n";
    for (auto pc : ParticleCombinationSet_) {
        std::cout << "  " << std::string(*pc);
        if (pc->parent()) {
            std::cout << "   \t in decay chain ";
            const ParticleCombination* parent = pc->parent();
            while (true) {
                if (parent->parent())
                    parent = parent->parent();
                else
                    break;
            }
            std::cout << std::string(*parent);
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
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
    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // check parent
    if (! ParticleCombination::equivUpButLambda(A->sharedParent(), B->sharedParent()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpAndDownButLambda::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivDownButLambda(A, B))
        return false;

    // check parent
    if (! ParticleCombination::equivUpButLambda(A->sharedParent(), B->sharedParent()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpAndDown::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
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
    // compare shared_ptr addresses
    // if both are nullptr, also return true
    if (A == B)
        return true;

    // if one is null_ptr return false
    if (A == nullptr or B == nullptr)
        return false;

    if (!ParticleCombination::equivByOrderlessContent(A->parent(), B->parent()))
        return false;

    return operator()(A->parent(), B->parent());
}

}
