#include "StatusManager.h"

#include "CalculationStatus.h"
#include "Parameter.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
StatusManager::StatusManager(const DataAccessorSet& sDA)
    : Statuses_(sDA.size())
{
    for (const auto& da : sDA) {
        Statuses_[da->index()].resize(da->CachedValues().size());
        for (const auto& cdv : da->CachedValues())
            Statuses_[da->index()][cdv->index()].resize(da->nSymmetrizationIndices());
    }
}

//-------------------------
CachedValue::Status& StatusManager::status(const CachedValue& cdv, size_t sym_index)
{
    return status(cdv.owner()->index(), cdv.index(), sym_index);
}

}
