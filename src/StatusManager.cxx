#include "StatusManager.h"

#include "CachedDataValue.h"
#include "DataAccessor.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
StatusManager::StatusManager(const DataAccessorSet& sDA)
    : Statuses_(sDA.size())
{
    for (const auto& da : sDA) {
        Statuses_[da->index()].resize(da->cachedDataValues().size());
        for (const auto& cdv : da->cachedDataValues())
            Statuses_[da->index()][cdv->index()].resize(da->maxSymmetrizationIndex() + 1);
    }
}

//-------------------------
CachedDataValue::Status& StatusManager::status(const CachedDataValue& cdv, size_t sym_index)
{
    return status(cdv.owner()->index(), cdv.index(), sym_index);
}

// //-------------------------
// void StatusManager::set(std::shared_ptr<CachedDataValue> cdv, const VariableStatus& stat)
// {
//     for (auto& s : Statuses_[cdv->owner()->index()][cdv->index()])
//         s.Variable = stat;
// }

//-------------------------
void StatusManager::set(const CachedDataValue& cdv, const CalculationStatus& stat)
{
    for (auto& s : Statuses_[cdv.owner()->index()][cdv.index()])
        s.Calculation = stat;
}

//-------------------------
void StatusManager::set(const DataAccessor& da, const CalculationStatus& stat)
{
    for (auto& v : Statuses_[da.index()])
        for (auto& s : v)
            s.Calculation = stat;
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

//-------------------------
void StatusManager::setAll(const VariableStatus& stat)
{
    for (auto& v1 : Statuses_)
        for (auto& v2 : v1)
            for (auto& s : v2)
                s.Variable = stat;
}

//-------------------------
void StatusManager::updateCalculationStatus(const CachedDataValue& cdv, const std::shared_ptr<ParticleCombination>& pc, size_t sym_index)
{
    // if already marked uncalculated, continue
    if (status(cdv, sym_index) == CalculationStatus::uncalculated)
        return;

    // check CachedDataValue dependencies
    for (const auto& c : cdv.cachedDataValueDependencies()) {

        if (c->owner() == cdv.owner()) {

            updateCalculationStatus(*c, pc, sym_index);

            if (status(*c, sym_index) == CalculationStatus::uncalculated or status(*c, sym_index) == VariableStatus::changed) {
                status(cdv, sym_index) = CalculationStatus::uncalculated;
                return;
            }
        } else if (c->owner()->hasParticleCombination(pc)) {
            auto s = c->owner()->symmetrizationIndex(pc);
            updateCalculationStatus(*c, pc, s);
            if (status(*c, s) == CalculationStatus::uncalculated or status(*c, s) == VariableStatus::changed) {
                status(cdv, sym_index) = CalculationStatus::uncalculated;
                return;
            }
        }
    }

    // if changed above, continue
    if (status(cdv, sym_index) == CalculationStatus::uncalculated)
        return;

    // check CachedDataValue daughter dependencies
    for (const auto dcdv : cdv.daughterCachedDataValueDependencies()) {

        const auto& dpc = pc->daughters()[dcdv.Daughter];

        if (dcdv.CDV->owner() != cdv.owner() and !dcdv.CDV->owner()->hasParticleCombination(dpc))
            continue;

        auto s = dcdv.CDV->owner()->symmetrizationIndex(dpc);

        updateCalculationStatus(*dcdv.CDV, dpc, s);

        if (status(*dcdv.CDV, s) == CalculationStatus::uncalculated or status(*dcdv.CDV, s) == VariableStatus::changed) {
            status(cdv, sym_index) = CalculationStatus::uncalculated;
            return;
        }
    }
}

//-------------------------
void StatusManager::updateCalculationStatuses(const CachedDataValue& cdv)
{
    // check parameter dependencies first,
    // since they don't depend on ParticleCombinations
    for (const auto& p : cdv.parameterDependencies())
        if (p->variableStatus() == VariableStatus::changed) {
            set(cdv, CalculationStatus::uncalculated);
            return;
        }

    /////////////////////////
    // check ParticleCombination-dependent dependencies
    for (const auto& kv : cdv.owner()->symmetrizationIndices())
        updateCalculationStatus(cdv, kv.first, kv.second);
}

//-------------------------
void StatusManager::updateCalculationStatuses(const DataAccessorSet& sDA)
{
    for (const auto& da : sDA)
        for (const auto& cdv : da->cachedDataValues())
            updateCalculationStatuses(*cdv);
}

}
