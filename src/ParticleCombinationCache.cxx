#include "ParticleCombinationCache.h"

#include "Exceptions.h"
#include "logging.h"

#include <algorithm>

namespace yap {

//-------------------------
/// \todo Shouldn't this management be done in the creation process?
ParticleCombinationCache::ParticleCombinationCache(std::vector<shared_ptr_type> V)
    : WeakPtrCache(V)
{
    // set parents from isp down
    for (auto& pc : V)
        setLineage(pc);

    // check consistency
    if (!consistent())
        throw exceptions::Exception("Inconsistent ParticleCombinationCache", "ParticleCombinationCache::ParticleCombinationCache");

    // remove expired cache entries
    removeExpired();
}

//-------------------------
void ParticleCombinationCache::setLineage(shared_ptr_type pc)
{
    /// \todo why do we need to do this here?
    for (auto& D : pc->Daughters_) {
        // copy D
        auto d = std::make_shared<type>(*D);
        // set copy's parent to pc
        std::const_pointer_cast<ParticleCombination>(d)->Parent_ = pc;
        // call recursively
        setLineage(pc);
        // replace daughter with copy (from cache)
        std::const_pointer_cast<ParticleCombination>(D) = std::const_pointer_cast<ParticleCombination>(operator[](d));
        // std::const_pointer_cast<ParticleCombination>(D).swap(std::const_pointer_cast<ParticleCombination>(operator[](d)));
    }
}

//-------------------------
ParticleCombinationCache::weak_ptr_type ParticleCombinationCache::find(const std::vector<ParticleIndex>& I) const
{
    std::vector<shared_ptr_type> V;
    V.reserve(I.size());
    for (auto& i : I) {
        auto a = find(std::make_shared<type>(i));
        if (a.expired())
            return weak_ptr_type();
        V.push_back(a.lock());
    }
    return find(std::make_shared<type>(V));
}

//-------------------------
void ParticleCombinationCache::addToCache(shared_ptr_type pc)
{
    // add daughters into cache, too (recursively)
    for (auto& d : pc->daughters())
        addToCache(d);
    // add self into cache
    WeakPtrCache::addToCache(pc);
    // setLineage(pc);
}

//-------------------------
ParticleCombinationCache::shared_ptr_type ParticleCombinationCache::operator[](const std::vector<ParticleIndex>& I)
{
    std::vector<shared_ptr_type> V;
    for (ParticleIndex i : I)
        V.push_back(operator[](i));
    return operator[](V);
}

//-------------------------
bool ParticleCombinationCache::consistent() const
{
    bool C = true;
    for (auto& wpc : *this) {
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
        s += to_string(*pc);
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
