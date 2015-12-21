#include "ParticleCombinationCache.h"

#include <cassert>

namespace yap {

//-------------------------
/// \todo Shouldn't this management be done in the creation process?
ParticleCombinationCache::ParticleCombinationCache(std::vector<std::shared_ptr<ParticleCombination> > ispPCs)
{
    // add particle combinations to cache
    for (auto& pc : ispPCs)
        operator[](pc);

    // set parents from isp down
    for (auto& pc : ispPCs)
        setLineage(pc);

    // check consistency
    assert(consistent());

    // remove expired cache entries
    removeExpired();
}

//-------------------------
void ParticleCombinationCache::setLineage(std::shared_ptr<ParticleCombination> pc)
{
    /// \todo why do we need to do this here?
    for (auto& D : pc->Daughters_) {
        // copy D
        auto d = std::make_shared<ParticleCombination>(*D);
        // set copy's parent to pc
        d->Parent_ = pc;
        // call recursively
        setLineage(pc);
        // replace daughter with copy (from cache)
        // D.swap(operator[](d));
        D = operator[](d);
    }
}

//-------------------------
ParticleCombinationCache::cache_type::key_type ParticleCombinationCache::find(const ParticleCombination* pc) const
{
    if (pc == nullptr)
        return cache_type::key_type();

    auto it = std::find_if(Cache_.begin(), Cache_.end(), [&](const cache_type::key_type & wpc) {return wpc.lock().get() == pc;});

    if (it == Cache_.end())
        return cache_type::key_type();

    return *it;
}

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombinationCache::operator[](std::shared_ptr<const ParticleCombination> pc)
{
    auto wpc = find(pc);

    // if wpc is valid
    if (!wpc.expired())
        return wpc.lock();

    // else add weak pointer to Cache_
    Cache_.emplace(pc);
    // and return original shared pointer
    return pc;
}

//-------------------------
std::shared_ptr<const ParticleCombination> ParticleCombinationCache::operator[](std::vector<ParticleIndex> I)
{
    ParticleCombinationVector V;
    for (ParticleIndex i : I)
        V.push_back(operator[](i));
    return operator[](V);
}

//-------------------------
void ParticleCombinationCache::removeExpired()
{
    for (auto it = Cache_.begin(); it != Cache_.end(); ) {
        if (it->expired())
            it = Cache_.erase(it);
        else
            it++;
    }
}

//-------------------------
bool ParticleCombinationCache::consistent() const
{
    bool C = true;
    for (auto& wpc : Cache_) {
        if (wpc.expired())
            continue;
        C &= wpc.lock()->consistent();
    }
    return C;
}

//-------------------------
std::string to_string(const ParticleCombinationCache& C)
{
    std::string s;
    for (auto wpc : C) {
        auto pc = wpc.lock();
        s += std::string(*pc);
        auto pt = pc->parent();
        if (pt) {
            s += " in decay chain ";
            while (pt->parent())
                pt = pt->parent();
            s += to_string(*pt);
        }
        s += "\n";
    }
    s.erase(s.size() - 1, 1);
    return s;
}

}
