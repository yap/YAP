#include "DataAccessor.h"

#include "DataPoint.h"
#include "Model.h"
#include "logging.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(ParticleCombination::Equiv* equiv) :
    ReportsModel(),
    ReportsParticleCombinations(),
    Equiv_(equiv),
    Size_(0),
    Index_(-1)
{
}

//-------------------------
int DataAccessor::maxSymmetrizationIndex() const
{
    /// I don't know why, but std::max_element returns wrong numbers sometimes!
    //return std::max_element(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), SymmetrizationIndices_.value_comp())->second;
    int max(-1);
    for (auto& kv : SymmetrizationIndices_)
        if (int(kv.second) > max)
            max = kv.second;

    return max;
}

//-------------------------
ParticleCombinationVector DataAccessor::particleCombinations() const
{
    ParticleCombinationVector particleCombinations;
    particleCombinations.reserve(SymmetrizationIndices_.size());

    for (auto& kv : SymmetrizationIndices_)
        particleCombinations.push_back(kv.first);

    return particleCombinations;
}

//-------------------------
void DataAccessor::printParticleCombinations() const
{
    LOG(INFO) << data_accessor_type();
    for (auto& kv : SymmetrizationIndices_) {
        auto p = kv.first;
        while (p->parent())
            p = p->parent();
        if (p == kv.first)
            LOG(INFO) << kv.second << " : " << *(kv.first);
        else
            LOG(INFO) << kv.second << " : " << *(kv.first) << ". In " << *p;
    }
}

//-------------------------
bool DataAccessor::hasParticleCombination(const std::shared_ptr<ParticleCombination>& c,
        const ParticleCombination::Equiv& equiv) const
{
    auto it = std::find_if(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(),
                           [&](const ParticleCombinationMap<unsigned>::value_type & kv)
    { return equiv(kv.first, c); } );
    return it != SymmetrizationIndices_.end();
}

//-------------------------
bool DataAccessor::consistent() const
{
    if (SymmetrizationIndices_.empty()) {
        FLOG(ERROR) << "SymmetrizationIndices_ is empty.";
        return false;
    }

    bool result = true;

    for (auto& kv : SymmetrizationIndices_)
        result &= kv.first->consistent();

    // check CachedDataValues_
    for (auto& c : CachedDataValues_) {
        if (c->CalculationStatus_.size() == 0) {
            FLOG(ERROR) << "c->CalculationStatus_.size() == 0";
            result = false;
        }
        if (c->CalculationStatus_.size() != c->VariableStatus_.size()) {
            FLOG(ERROR) << "c->CalculationStatus_.size() != c->VariableStatus_.size()";
            result = false;
        }

        for (unsigned i = 0; i < c->CalculationStatus_.size(); ++i) {
            if (int(c->CalculationStatus_[i].size()) != maxSymmetrizationIndex() + 1) {
                FLOG(ERROR) << "c->CalculationStatus_[i].size() != maxSymmetrizationIndex() + 1 ("
                            << c->CalculationStatus_.size() << " != " << maxSymmetrizationIndex() + 1 << ")";
                DEBUG("c's Owner " << c->owner() << " " << typeid(*c->owner()).name());
                result = false;
            }
            if (int(c->VariableStatus_[i].size()) != maxSymmetrizationIndex() + 1) {
                LOG(ERROR) << "DataAccessor::consistent() - c->VariableStatus_[i].size() != maxSymmetrizationIndex() + 1";
                LOG(ERROR) << c->VariableStatus_.size() << " != " << maxSymmetrizationIndex() + 1;
                DEBUG("c's Owner " << c->owner() << " " << typeid(*c->owner()).name());
                result = false;
            }
        }

    }

    return result;
}

//-------------------------
unsigned DataAccessor::addParticleCombination(std::shared_ptr<ParticleCombination> c)
{
    if (!c)
        throw exceptions::Exception("ParticleCombination empty", "DataAccessor::addParticleCombination");

    if (hasParticleCombination(c))
        // c is already in map
        return symmetrizationIndex(c);

    // check to see if new member equates to existing member
    for (auto& kv : SymmetrizationIndices_)
        if ((*Equiv_)(kv.first, c)) {
            // equating member found; set index; return
            SymmetrizationIndices_[c] = kv.second;
            return kv.second;
        }

    // else assign to current size = highest current index + 1
    unsigned size = maxSymmetrizationIndex() + 1;
    SymmetrizationIndices_[c] = size;

    for (auto& c : CachedDataValues_)
        c->setNumberOfSymmetrizations(size + 1);

    return size;
}

//-------------------------
void DataAccessor::pruneSymmetrizationIndices()
{
    if (!model())
        throw exceptions::Exception("Model not set", "DataAccessor::pruneSymmetrizationIndices");

    // remove entries that don't trace back to the ISP
    for (auto it = SymmetrizationIndices_.begin(); it != SymmetrizationIndices_.end(); ) {
        // find the top-most parent
        auto pc = it->first;
        while (pc->parent())
            pc = pc->parent();
        // check if it's not an ISP
        if (pc->indices().size() != model()->finalStateParticles().size())
            // erase
            it = SymmetrizationIndices_.erase(it);
        else
            it++;
    }

    // fix indices now for holes

    // collect used indices
    std::set<unsigned> used;
    for (const auto& kv : SymmetrizationIndices_)
        used.insert(kv.second);

    // repair
    unsigned index = 0;
    while (index < used.size()) {

        // if index is not used
        if (used.find(index) == used.end()) {
            // clear used
            used.clear();
            // reduce all indices greater than index by 1
            // and rebuild used
            for (auto& kv : SymmetrizationIndices_) {
                if (kv.second > index)
                    kv.second -= 1;
                used.insert(kv.second);
            }
        }

        //if index is (now) used, increment by 1
        if (used.find(index) != used.end())
            index += 1;

    }

    unsigned size = maxSymmetrizationIndex() + 1;
    for (auto& c : CachedDataValues_)
        c->setNumberOfSymmetrizations(size);
}

//-------------------------
void DataAccessor::addToModel()
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataAccessor::addToModel");
    model()->addDataAccessor(this);
}

//-------------------------
void DataAccessor::updateGlobalCalculationStatuses()
{
    for (auto& c : CachedDataValues_) {
        for (auto& kv : SymmetrizationIndices_) {
            DEBUG("updateGlobalCalculationStatuses for " << typeid(*this).name() << " " << dynamic_cast<DataAccessor*>(this) << " for " << * (kv.first));
            c->updateGlobalCalculationStatus(kv.first, kv.second);
        }
    }
}

//-------------------------
void DataAccessor::setNumberOfDataPartitions(unsigned n)
{
    for (auto& c : CachedDataValues_)
        c->setNumberOfDataPartitions(n);
}

//-------------------------
void DataAccessor::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    for (auto& c : CachedDataValues_)
        c->resetCalculationStatus(dataPartitionIndex);
}

//-------------------------
void DataAccessor::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    for (auto& c : CachedDataValues_)
        c->setVariableStatus(kUnchanged, dataPartitionIndex);
}

//-------------------------
void DataAccessor::setParameterFlagsToUnchanged()
{
    for (auto& c : CachedDataValues_)
        for (auto& p : c->ParametersItDependsOn_)
            if (p->variableStatus() == kChanged)
                p->setVariableStatus(kUnchanged);
}

//-------------------------
void removeExpired(DataAccessorSet& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}


}

