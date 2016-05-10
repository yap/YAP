#include "DataAccessor.h"

#include "CachedDataValue.h"
#include "logging.h"
#include "Model.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(ParticleCombination::Equiv* equiv) :
    ReportsParticleCombinations(),
    Equiv_(equiv),
    NIndices_(0),
    Size_(0),
    Index_(-1)
{
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

    // object for recording successing of emplacement
    auto it_b = std::make_pair(SymmetrizationIndices_.end(), false);

    // check to see if new member equates to existing member
    for (auto& kv : SymmetrizationIndices_)
        if ((*Equiv_)(kv.first, c))
            // equating member found; set index; return
            it_b = SymmetrizationIndices_.emplace(c, kv.second);
    // if c is new but equates to existing member, it_b.first != end and it_b.second is true
    // if c is existing member, it_b.first != end and it_b.second is false
    // if c does not equate to existing member, it_b.first = end, it_b.second is false

    // else assign to one higher than current highest index
    if (it_b.first == SymmetrizationIndices_.end()) {
        it_b = SymmetrizationIndices_.emplace(c, static_cast<unsigned>(NIndices_));
        // if emplacement successful, it_b.first != end, it_b.second = true

        // if successfully emplaced, increase NIndices_
        if (it_b.second)
            ++NIndices_;
    }

    if (it_b.first == SymmetrizationIndices_.end())
        throw exceptions::Exception("Failed to emplace new SymmetrizationIndices element.", "DataAccessor::addParticleCombination");

    if (it_b.second)
        // add to ParticleCombinations_
        ParticleCombinations_.push_back(c);

    // return c's index in map
    return it_b.first->second;
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

    // reset NIndices_
    NIndices_ = 0;
    for (const auto& kv : SymmetrizationIndices_)
        NIndices_ = std::max(kv.second + 1, NIndices_);
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

