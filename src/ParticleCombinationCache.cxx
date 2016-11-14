#include "ParticleCombinationCache.h"

#include "container_utils.h"
#include "Exceptions.h"

#include <algorithm>
#include <cmath>

namespace yap {

//-------------------------
ParticleCombinationCache::shared_ptr_type ParticleCombinationCache::create_copy(const ParticleCombination& other) const
{
    if (is_final_state_particle_combination(other))
        return create_fsp(other.indices()[0]);

    auto pc = shared_ptr_type(new ParticleCombination());

    // copy daughters (recursively) and add them
    for (auto& d : other.Daughters_)
        std::const_pointer_cast<non_const_type>(pc)->addDaughter(*std::const_pointer_cast<non_const_type>(create_copy(*d)));

    return pc;
}

//-------------------------
ParticleCombinationCache::shared_ptr_type ParticleCombinationCache::create_composite(const ParticleCombinationVector& D) const
{
    auto pc = shared_ptr_type(new ParticleCombination());

    for (auto& d : D)
        std::const_pointer_cast<non_const_type>(pc)->addDaughter(*std::const_pointer_cast<non_const_type>(create_copy(*d)));

    return pc;
}

//-------------------------
void ParticleCombinationCache::addToCache(shared_ptr_type pc)
{
    // add daughters into cache, too (recursively)
    for (auto& d : pc->daughters()) {
        // look for daughter in cache
        auto w = find(d);
        // if it's not there, add it
        if (w.expired())
            addToCache(d);
    }
    // add self into cache
    WeakPtrCache::addToCache(pc);
    // setLineage(pc);
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

}
