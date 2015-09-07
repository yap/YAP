#include "DataAccessor.h"

#include "DataPoint.h"
#include "logging.h"

namespace yap {

unsigned DataAccessor::GlobalIndex = 0;

//-------------------------
DataAccessor::DataAccessor(ParticleCombination::Equiv equiv) :
    Equiv_(equiv),
    Index_(0)
{
    // assign a running index to this DataAccessor
    /// \todo Come up with something smarter
    Index_ = GlobalIndex++;
}

//-------------------------
DataAccessor::DataAccessor(const DataAccessor& other) :
    Equiv_(other.Equiv_),
    CalculationStatuses_(other.CalculationStatuses_),
    SymmetrizationIndices_(other.SymmetrizationIndices_),
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

    // find number of indices:
    unsigned n_indices = std::max_element(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), SymmetrizationIndices_.value_comp())->second;
    // check that CalculationStatuses_ has right number of entries
    if (CalculationStatuses_.size() != n_indices + 1) {
        LOG(ERROR) << "DataAccessor::consistent() - CalculationStatuses_ has wrong number of entries ("
                   << CalculationStatuses_.size() << " != " << n_indices + 1 << ")";
        result = false;
    }

    for (auto& kv : SymmetrizationIndices_)
        result &= kv.first->consistent();

    // find number of indices:

    return result;
}

//-------------------------
void DataAccessor::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
    if (SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end())
        // c is already in map
        return;

    // check to see if new member equates to existing member
    for (auto& kv : SymmetrizationIndices_)
        if (Equiv_(kv.first, c)) {
            // equating member found; set index; return
            SymmetrizationIndices_[c] = kv.second;
            return;
        }

    // else assign to current size = highest current index + 1
    SymmetrizationIndices_[c] = CalculationStatuses_.size();
    // and add entry to
    CalculationStatuses_.push_back(kUncalculated);

}

}

