#include "ParticleCombination.h"

#include "logging.h"

#include <algorithm>

namespace yap {

//-------------------------
ParticleCombination::ParticleCombination(ParticleIndex index) :
    Indices_(1, index)
{
}

//-------------------------
bool ParticleCombination::addDaughter(std::shared_ptr<ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        LOG(ERROR) << "ParticleCombination::addDaughter - daughter contains no indices.";
        return false;
    }

    Daughters_.push_back(daughter);
    Indices_.insert(Indices_.end(), daughter->indices().begin(), daughter->indices().end());
    return true;
}

//-------------------------
bool ParticleCombination::consistent() const
{
    bool result = true;

    // create unique_copy of Indices_
    std::vector<ParticleIndex> U;
    std::unique_copy(Indices_.begin(), Indices_.end(), U.begin());
    // check unique_copy-object's size == this object's size
    if (U.size() != Indices_.size()) {
        LOG(ERROR) << "ParticleCombination::consistent - index vector contains duplicate entries.";
        result = false;
    }

    // check daughers
    for (unsigned i = 0; i < Daughters_.size(); ++i)
        result &= Daughters_[i]->consistent();

    return result;
}

//-------------------------
bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{
    // Check indices
    if (A.Indices_.size() != B.Indices_.size())
        return false;
    if (!std::equal(A.Indices_.begin(), A.Indices_.end(), B.Indices_.begin()))
        return false;

    // Check daughters
    if (A.Daughters_.size() != B.Daughters_.size())
        return false;
    for (unsigned i = 0; i < A.Daughters_.size(); ++i)
        if (*(A.Daughters_[i]) != *(B.Daughters_[i]))
            return false;

    // a match!
    return true;
}

}
