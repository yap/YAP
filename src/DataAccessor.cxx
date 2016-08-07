#include "DataAccessor.h"

#include "CachedValue.h"
#include "logging.h"
#include "Model.h"
#include "RequiresHelicityAngles.h"
#include "RequiresMeasuredBreakupMomenta.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(const ParticleCombinationEqualTo& equal) :
    Equal_(equal),
    NIndices_(0),
    Size_(0),
    Index_(-1)
{
}

//-------------------------
void DataAccessor::printParticleCombinations() const
{
    for (auto& kv : SymmetrizationIndices_)
        LOG(INFO) << kv.second << " : " << to_string_with_parent(*(kv.first));
}

//-------------------------
bool DataAccessor::consistent() const
{
    bool C = true;

    // check CachedValues_
    for (auto& c : CachedValues_)
        if (c->owner() != this) {
            FLOG(ERROR) << "CachedValue's owner != this";
            C &= false;
        }

    return C;
}

//-------------------------
void DataAccessor::addParticleCombination(std::shared_ptr<ParticleCombination> c)
{
    if (!c)
        throw exceptions::Exception("ParticleCombination empty", "DataAccessor::addParticleCombination");

    // check if c is already key in SymmetrizationIndices_
    if (SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end())
        return;

    // search for match using Equal
    auto it = std::find_if(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(),
                           [&](const ParticleCombinationMap<unsigned>::value_type & kv)
    {return Equal_(c, kv.first);});

    // if found, use found index
    if (it != SymmetrizationIndices_.end()) {
        SymmetrizationIndices_.emplace(c, it->second);
        return;
    }

    // else assign to one higher than current highest index
    SymmetrizationIndices_.emplace(c, static_cast<unsigned>(NIndices_));
    // and increase current highest index
    ++NIndices_;
}

//-------------------------
void DataAccessor::pruneSymmetrizationIndices()
{
    if (!model())
        throw exceptions::Exception("Model not set", "DataAccessor::pruneSymmetrizationIndices");

    // remove entries that don't trace back to an ISP
    for (auto it = SymmetrizationIndices_.begin(); it != SymmetrizationIndices_.end(); ) {
        if (is_from_initial_state_particle_combination(*it->first, *model()))
            ++it;
        else
            it = SymmetrizationIndices_.erase(it);
    }

    //
    // fix indices now for holes
    //

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

    // reset NIndices_
    NIndices_ = 0;
    for (const auto& kv : SymmetrizationIndices_)
        NIndices_ = std::max(kv.second + 1, NIndices_);
}

//-------------------------
void DataAccessor::registerWithModel()
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataAccessor::registerWithModel");
    
    if (model()->locked())
        throw exceptions::Exception("Model is locked and cannot be modified", "DataAccessor::registerWithModel");
    
    // if HelicityAngles is required
    if (dynamic_cast<RequiresHelicityAngles*>(this) and dynamic_cast<RequiresHelicityAngles*>(this)->requiresHelicityAngles())
        const_cast<Model*>(model())->requireHelicityAngles();
    
    // if MeasuredBreakupMomenta is required
    if (dynamic_cast<RequiresMeasuredBreakupMomenta*>(this) and dynamic_cast<RequiresMeasuredBreakupMomenta*>(this)->requiresMeasuredBreakupMomenta())
        const_cast<Model*>(model())->requireMeasuredBreakupMomenta();
    
    // if stores nothing, do nothing
    if (size() == 0)
        return;

    // insert this into Model's DataAccessors_
    const_cast<Model*>(model())->DataAccessors_.insert(this);
}

//-------------------------
void DataAccessor::addCachedValue(std::shared_ptr<CachedValue> c)
{
    if (not c)
        throw exceptions::Exception("CachedValue is NULL", "DataAccessor::addCachedValue");

    // add CachedValue
    if (CachedValues_.insert(c).second) {
        // if insertion was successful

        // set its index
        c->setIndex(CachedValues_.size() - 1);

        // set its position
        c->setPosition(size());

        // increase data size to accommodate CachedValue
        increaseSize(c->size());
    }
}

//-------------------------
void remove_expired(DataAccessorSet& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}


}

