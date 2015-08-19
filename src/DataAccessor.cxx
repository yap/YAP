#include "DataAccessor.h"

#include "logging.h"

namespace yap {

unsigned DataAccessor::GlobalIndex = 0;

//-------------------------
DataAccessor::DataAccessor() :
    Recalculate_(true),
    Index_(0)
{
    // assign a running index to this DataAccessor
    /// \todo Come up with something smarter
    Index_ = GlobalIndex++;
}

//-------------------------
DataAccessor::DataAccessor(const DataAccessor& other) :
    Recalculate_(other.Recalculate_),
    Index_(0)
{
    // uses new index
    Index_ = GlobalIndex++;
}

//-------------------------
std::vector<std::shared_ptr<ParticleCombination> > DataAccessor::particleCombinations() const
{
    std::vector<std::shared_ptr<ParticleCombination> > retVal;
    for (auto& kv : SymmetrizationIndices_)
        retVal.push_back(kv.first);

    return retVal;
}

//-------------------------
bool DataAccessor::consistent() const
{
    if (SymmetrizationIndices_.empty()) {
        LOG(ERROR) << "DataAccessor::consistent() - SymmetrizationIndices_ is empty.";
        return false;
    }

    bool result = true;
    for (auto& kv : SymmetrizationIndices_)
        result &= kv.first->consistent();

    return result;
}

//-------------------------
void DataAccessor::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
    if (SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end())
        // c is already in map
        return;

    // if map empty, add
    if (SymmetrizationIndices_.empty()) {
        SymmetrizationIndices_[c] = 0;
        return;
    }

    // else check to see if new member equates to existing member
    // and search for highest index otherwise
    unsigned index = 0;
    for (auto& kv : SymmetrizationIndices_) {
        if (areEqual(kv.first, c)) {
            // equating member found; set index; return
            SymmetrizationIndices_[c] = kv.second;
            return;
        }
        index = std::max(index, kv.second);
    }

    SymmetrizationIndices_[c] = (index + 1);
}

}

