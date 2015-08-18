#include "DataAccessor.h"

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
void DataAccessor::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
  if (SymmetrizationIndices_.find(c) == SymmetrizationIndices_.end()) {
    // simple running index
    unsigned index = SymmetrizationIndices_.size();
    SymmetrizationIndices_[c] = index;
  }
}

}

