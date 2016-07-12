#include "StatusManager.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataAccessor.h"
#include "Exceptions.h"
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

//-------------------------
void StatusManager::copyCalculationStatuses(const StatusManager& sm)
{
    if (Statuses_.size() != sm.Statuses_.size())
        throw exceptions::Exception("size mismatch", "StatusManager::setCalculationStatus");
    for (size_t i = 0; i < Statuses_.size(); ++i) {
        if (Statuses_[i].size() != sm.Statuses_[i].size())
            throw exceptions::Exception("size mismatch", "StatusManager::setCalculationStatus");
        for (size_t j = 0; j < Statuses_[i].size(); ++j) {
            if (Statuses_[i][j].size() != sm.Statuses_[i][j].size())
                throw exceptions::Exception("size mismatch", "StatusManager::setCalculationStatus");
            for (size_t k = 0; k < Statuses_[i][j].size(); ++k)
                Statuses_[i][j][k].Calculation = sm.Statuses_[i][j][k].Calculation;
        }
    }
}

}
