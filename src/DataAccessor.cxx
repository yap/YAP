#include "DataAccessor.h"

#include "DataPoint.h"
#include "InitialStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(ParticleCombination::Equiv* equiv) :
    ReportsInitialStateParticle(),
    ReportsParticleCombinations(),
    Equiv_(equiv),
    Size_(0),
    Index_(0)
{
    // index is set later by InitialStateParticle::setDataAcessorIndices()
    // via InitialStateParticle::prepare()
}

//-------------------------
DataAccessor::~DataAccessor()
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
        LOG(ERROR) << "DataAccessor::consistent() - SymmetrizationIndices_ is empty.";
        return false;
    }

    bool result = true;

    for (auto& kv : SymmetrizationIndices_)
        result &= kv.first->consistent();

    // check CachedDataValues_
    for (CachedDataValue* c : CachedDataValues_) {
        if (c->CalculationStatus_.size() == 0) {
            LOG(ERROR) << "DataAccessor::consistent() - c->CalculationStatus_.size() == 0";
            result = false;
        }
        if (c->CalculationStatus_.size() != c->VariableStatus_.size()) {
            LOG(ERROR) << "DataAccessor::consistent() - c->CalculationStatus_.size() != c->VariableStatus_.size()";
            result = false;
        }

        for (unsigned i = 0; i < c->CalculationStatus_.size(); ++i) {
            if (int(c->CalculationStatus_[i].size()) != maxSymmetrizationIndex() + 1) {
                LOG(ERROR) << "DataAccessor::consistent() - c->CalculationStatus_[i].size() != maxSymmetrizationIndex() + 1";
                LOG(ERROR) << c->CalculationStatus_.size() << " != " << maxSymmetrizationIndex() + 1;
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
void DataAccessor::addParticleCombination(std::shared_ptr<ParticleCombination> c)
{
    if (hasParticleCombination(c))
        // c is already in map
        return;

    // check to see if new member equates to existing member
    for (auto& kv : SymmetrizationIndices_)
        if ((*Equiv_)(kv.first, c)) {
            // equating member found; set index; return
            SymmetrizationIndices_[c] = kv.second;
            return;
        }

    // else assign to current size = highest current index + 1
    unsigned size = maxSymmetrizationIndex() + 1;
    SymmetrizationIndices_[c] = size;

    for (CachedDataValue* d : CachedDataValues_)
        d->setNumberOfSymmetrizations(size + 1);
}

//-------------------------
void DataAccessor::addCachedDataValue(CachedDataValue* c)
{
    c->setNumberOfSymmetrizations(maxSymmetrizationIndex() + 1);
    CachedDataValues_.insert(c);
}

//-------------------------
std::vector<double>& DataAccessor::data(DataPoint& d, unsigned i) const
{
    // dynamically allocate memory as needed
    if (d.Data_.size() <= Index_)
        d.Data_.resize(Index_ + 1);

    if (d.Data_[Index_].size() <= i)
        d.Data_[Index_].resize(i + 1);

    return d.Data_[Index_][i];
}

//-------------------------
const std::vector<double>& DataAccessor::data(const DataPoint& d, unsigned i) const
{
#ifdef ELPP_DISABLE_DEBUG_LOGS
    return d.Data_[Index_][i];
#else
    return d.Data_.at(Index_).at(i);
#endif
}


//-------------------------
void DataAccessor::updateGlobalCalculationStatuses()
{
    for (CachedDataValue* c : CachedDataValues_) {
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

}

