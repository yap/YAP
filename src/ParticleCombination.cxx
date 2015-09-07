#include "ParticleCombination.h"

#include "logging.h"

#include <algorithm>
#include <assert.h>

namespace yap {

//-------------------------
ParticleCombination::ParticleCombination()
{
}

//-------------------------
ParticleCombination::ParticleCombination(ParticleIndex index) :
    Indices_(1, index)
{
}

//-------------------------
ParticleCombination::ParticleCombination(std::vector<std::shared_ptr<ParticleCombination> > c)
{
    for (std::shared_ptr<ParticleCombination>& d : c)
        addDaughter(d);
}

//-------------------------
bool ParticleCombination::addDaughter(std::shared_ptr<ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        LOG(ERROR) << "ParticleCombination::addDaughter - daughter contains no indices.";
        return false;
    }

    /// \todo Check that new daughter does not share content with other daughters?



    /*if (daughter->Parent_ == nullptr) {
        // daughter has no parent yet. Set this as parent and add
        daughter->Parent_ = this;
    } else if (daughter->Parent_ == this) {
        // fine
    } else {
        LOG(ERROR) << "this should not happen!";

        // daughter has already different parent -> make copy, set this as parent and get unique shared_ptr
        std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*daughter));
        copy->Parent_ = this;
        std::shared_ptr<ParticleCombination> uniqueCopy = uniqueSharedPtr(copy);
        // need to swap so that parent of argument daughter is now set

        // \todo does not work since it would have to operate on the original shared_ptr object
        daughter.swap(uniqueCopy);
    }

    assert(daughter->Parent_ == this);*/

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

    // check daughters
    for (std::shared_ptr<ParticleCombination> d : Daughters_) {
        if (not d->parents().empty() and  std::find(d->parents().begin(), d->parents().end(), this) == d->parents().end()) {
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
bool ParticleCombination::sharesIndices(std::shared_ptr<ParticleCombination> B)
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
        /*if (daughter->Parent_ == nullptr) {
            // daughter has no parent yet. Set this as parent and add
            daughter->Parent_ = this;
        } else if (daughter->Parent_ == this) {
            // fine
        } else {
            // daughter has already different parent -> make copy, set this as parent and get unique shared_ptr
            std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*daughter));
            copy->Parent_ = this;
            std::shared_ptr<ParticleCombination> uniqueCopy = uniqueSharedPtr(copy);
            // need to swap so that parent of daughter is now set
            daughter.swap(uniqueCopy);
        }

        assert(daughter->Parent_ == this);*/

        daughter->addParent(this);
    }
}

//-------------------------
void ParticleCombination::addParent(ParticleCombination* parent)
{
    unsigned n = std::count(Parents_.begin(), Parents_.end(), parent);
    if (n == 0) {
        Parents_.push_back(parent);
        return;
    }

    if (n > 1) {
        LOG(ERROR) << "duplicate parent.";
    }
}

//-------------------------
bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{
    return ParticleCombination::equivUpAndDown(std::make_shared<ParticleCombination>(A), std::make_shared<ParticleCombination>(B));
}

/////////////////////////
// Static stuff:

std::set<std::shared_ptr<ParticleCombination> > ParticleCombination::ParticleCombinationSet_;

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(std::shared_ptr<ParticleCombination> pc)
{
    for (std::shared_ptr<ParticleCombination> d : ParticleCombinationSet_)
        if ((*pc) == (*d))
            return d;
    ParticleCombinationSet_.insert(pc);
    return pc;
}

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(ParticleIndex i)
{
    return uniqueSharedPtr(std::make_shared<ParticleCombination>(i));
}

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(std::vector<std::shared_ptr<ParticleCombination> > c)
{
    return uniqueSharedPtr(std::make_shared<yap::ParticleCombination>(c));
}

//-------------------------
void ParticleCombination::printParticleCombinationSet()
{
    std::cout << "ParticleCombination set:\n";
    for (auto pc : ParticleCombinationSet_) {
        std::cout << "  " << std::string(*pc) << "\n";
    }
    std::cout << std::endl;
}

//-------------------------
// Comparison shtuff:

ParticleCombination::Equiv ParticleCombination::equivBySharedPointer;
ParticleCombination::EquivDown ParticleCombination::equivDown;
ParticleCombination::EquivUpAndDown ParticleCombination::equivUpAndDown;
ParticleCombination::EquivByOrderedContent ParticleCombination::equivByOrderedContent;
ParticleCombination::EquivByOrderlessContent ParticleCombination::equivByOrderlessContent;

//-------------------------
bool ParticleCombination::EquivByOrderedContent::operator()(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B) const
{
    // Check indices
    if (A->indices().size() != B->indices().size())
        return false;
    if (!std::equal(A->indices().begin(), A->indices().end(), B->indices().begin()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivDown::operator()(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B) const
{
    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size())
        return false;
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        if (operator()(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUpAndDown::operator()(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B) const
{
    if (!ParticleCombination::equivDown(A, B))
        return false;

    /// \todo Can really only have one parent
    // check parents
    if (A->parents().size() != B->parents().size())
        return false;
    if (!std::equal(A->parents().begin(), A->parents().end(), B->parents().begin()))
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivByOrderlessContent::operator()(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B) const
{
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
