#include "DataAccessor.h"

#include "CachedDataValue.h"
#include "logging.h"
#include "Model.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(ParticleCombination::Equiv* equiv) :
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
    bool C = true;

    // check CachedDataValues_
    for (auto& c : CachedDataValues_)
        if (c->owner() != this) {
            FLOG(ERROR) << "CachedDataValue's owner != this";
            C &= false;
        }

    if (SymmetrizationIndices_.size() != ParticleCombinations_.size()) {
        FLOG(ERROR) << "SymmetrizationIndices_.size() " << SymmetrizationIndices_.size()
                    << " != ParticleCombinations_.size() " << ParticleCombinations_.size();
        C &= false;
    }

    return C;
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
            // and also add to ParticleCombinations_
            ParticleCombinations_.push_back(c);

            return kv.second;
        }

    // else assign to current size = highest current index + 1
    unsigned size = maxSymmetrizationIndex() + 1;
    SymmetrizationIndices_[c] = size;

    // and also add to ParticleCombinations_
    ParticleCombinations_.push_back(c);

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

    // (re)fill ParticleCombinations_
    ParticleCombinations_.clear();
    ParticleCombinations_.reserve(SymmetrizationIndices_.size());

    for (auto& kv : SymmetrizationIndices_)
        ParticleCombinations_.push_back(kv.first);
}

//-------------------------
void DataAccessor::addToModel()
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataAccessor::addToModel");
    const_cast<Model*>(static_cast<const DataAccessor*>(this)->model())->addDataAccessor(this);
}

//-------------------------
void DataAccessor::addCachedDataValue(std::shared_ptr<CachedDataValue> c)
{
    // add CachedDataValue
    if (CachedDataValues_.insert(c).second) {
        // if insertion was successful

        // set its index
        c->setIndex(CachedDataValues_.size() - 1);

        // set its position
        c->setPosition(size());

        // increase data size to accommodate CachedDataValue
        increaseSize(c->size());
    }
}

//-------------------------
void removeExpired(DataAccessorSet& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}

}

