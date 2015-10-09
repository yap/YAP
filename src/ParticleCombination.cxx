#include "ParticleCombination.h"

#include "logging.h"
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
ParticleCombination::ParticleCombination(std::vector<std::shared_ptr<const ParticleCombination> > c, char twoLambda) :
    Parent_(nullptr),
    TwoLambda_(twoLambda)
{
    for (auto& d : c)
        addDaughter(d);
}

//-------------------------
/*std::vector<std::shared_ptr<const ParticleCombination> > ParticleCombination::daughters() const
{
    std::vector<std::shared_ptr<const ParticleCombination> > daughters;
    for (auto& d : Daughters_)
        daughters.push_back(d);

    return daughters;
}*/

//-------------------------
const std::shared_ptr<const ParticleCombination> ParticleCombination::sharedParent() const
{
    if (! Parent_) {
        return std::shared_ptr<ParticleCombination>(Parent_);
    }

    for (auto& pc : ParticleCombinationSet_)
        if (Parent_ == pc.get())
            return pc;

    LOG(WARNING) << "ParticleCombination::parent() - could not find parent in ParticleCombinationSet_.";

    return std::shared_ptr<const ParticleCombination>(Parent_);
}

//-------------------------
bool ParticleCombination::addDaughter(std::shared_ptr<const ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        LOG(ERROR) << "ParticleCombination::addDaughter - daughter contains no indices.";
        return false;
    }

    /// Check that new daughter does not share content with other daughters?
    for (unsigned indexP : Indices_)
        for (unsigned indexD : daughter->indices())
            if (indexP == indexD) {
                LOG(ERROR) << "ParticleCombination::addDaughter - daughter contains indices that are already in parent.";
                return false;
            }

    // add daughter to vector
    Daughters_.push_back(daughter);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), daughter->indices().begin(), daughter->indices().end());

    return true;
}

//-------------------------
bool ParticleCombination::consistent() const
{
    // should have no daughters or 2 or more daughters:
    if (Daughters_.size() == 1) {
        LOG(ERROR) << "ParticleCombination::consistent() - has only one daughter.";
        return false;
    }

    // if empty, return inconsistent
    if (Indices_.empty()) {
        LOG(ERROR) << "ParticleCombination::consistent() - has no indices.";
        return false;
    }

    // if Daugthers_ empty, should have one and only index
    if (Daughters_.empty()) {
        if (Indices_.size() != 1) {
            LOG(ERROR) << "ParticleCombination::consistent() - contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
            return false;
        } else
            // don't need to check indices below
            return true;
    }

    // Check indices & and then daughters
    bool result = true;

    // create unique_copy of Indices_ (as set)
    std::set<ParticleIndex> U(Indices_.begin(), Indices_.end());
    // check unique_copy-object's size == this object's size
    if (U.size() != Indices_.size()) {
        LOG(ERROR) << "ParticleCombination::consistent - index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
        result = false;
    }

    // check if in ParticleCombinationSet_
    bool found(false);
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc.get() == this) {
            found = true;
            break;
        }
    }
    if (!found) {
        if (parent())
            LOG(ERROR) << "ParticleCombination::consistent - ParticleCombination is not in ParticleCombinationSet: " << std::string(*this) << " from decay " << std::string(*parent());
        else
            LOG(ERROR) << "ParticleCombination::consistent - ParticleCombination is not in ParticleCombinationSet: " << std::string(*this) << " (no parent)";
        result = false;
    }

    // check daughters
    for (auto& d : Daughters_) {
        if (d->parent() and  d->parent() != this) {
            LOG(ERROR) << "ParticleCombination::consistent - daughter's parent is not this ParticleCombination.";
            result = false;
        }
        result &= d->consistent();
    }

    return result;
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
void ParticleCombination::setParent(ParticleCombination* parent)
{
    Parent_ = parent;
}

//-------------------------
bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{
    return ParticleCombination::equivUpAndDown(std::make_shared<ParticleCombination>(A), std::make_shared<ParticleCombination>(B));
}

/////////////////////////
// Static stuff:

std::set<std::shared_ptr<const ParticleCombination> > ParticleCombination::ParticleCombinationSet_;

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
std::shared_ptr<const ParticleCombination> ParticleCombination::uniqueSharedPtr(std::vector<std::shared_ptr<const ParticleCombination> > c)
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

//-------------------------
bool ParticleCombination::EquivByOrderedContent::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivDownButLambda::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivDown::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivUpButLambda::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivUpAndDownButLambda::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivUpAndDown::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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
bool ParticleCombination::EquivByOrderlessContent::operator()(std::shared_ptr<const ParticleCombination> A, std::shared_ptr<const ParticleCombination> B) const
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

}
