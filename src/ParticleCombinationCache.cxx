#include "ParticleCombinationCache.h"

#include "container_utils.h"
#include "Exceptions.h"

#include <algorithm>
#include <cmath>
#include <iomanip>

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
ParticleCombinationCache::weak_ptr_type ParticleCombinationCache::findByUnorderedContent(const std::vector<unsigned>& I) const
{
    // look for entry with same content
    // by checking if indices() contains I and is the same size
    auto it = std::find_if(begin(), end(),
                           [&](const weak_ptr_type & w)
    {return !w.expired() and I.size() == w.lock()->indices().size() and contains(w.lock()->indices(), I);});
    if (it == end())
        // if not found
        return weak_ptr_type();

    return *it;
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

//-------------------------
std::ostream& print_daughters(std::ostream& os, const ParticleCombinationCache& C, const ParticleCombination& pc, unsigned ndig, std::string prefix, std::set<unsigned>& used)
{
    for (auto& d : pc.daughters()) {
        unsigned i = 0;
        for (auto& w : C) {
            if (!w.expired() and w.lock() == d) {
                os << std::setfill('0') << std::setw(ndig) << i << " : " << prefix << *w.lock();
                if (d->parent() != pc.shared_from_this())
                    os << " (parent unset)";
                os << std::endl;
                used.insert(i);
                break;
            }
            ++i;
        }
        print_daughters(os, C, *d, ndig, prefix + "    ", used);
    }
    return os;
}

//-------------------------
std::ostream& ParticleCombinationCache::print(std::ostream& os) const
{
    if (size() == 0)
        return os;

    // get number of digits for writing position in cache
    unsigned ndig = std::ceil(log10(size()));

    // find largest number of indices
    unsigned n_fsp = 0;
    for (auto& w : *this)
        if (!w.expired() and w.lock()->indices().size() > n_fsp)
            n_fsp = w.lock()->indices().size();

    std::set<unsigned> used;

    // print others as decay trees starting from ISP's
    unsigned i = 0;
    for (auto& w : *this) {
        if (!w.expired() and w.lock()->indices().size() == n_fsp) {
            // print isp
            os << std::setfill('0') << std::setw(ndig) << i << " : " << *w.lock() << std::endl;
            used.insert(i);

            // print daughters
            print_daughters(os, *this, *w.lock(), ndig, "    ", used);
        }

        ++i;
    }

    // print unused
    i = 0;
    for (auto& w : *this) {
        if (!w.expired() and used.find(i) == used.end())
            os << std::setfill('0') << std::setw(ndig) << i << " : (unused, " << w.use_count() << ") " << *w.lock() << std::endl;
        ++i;
    }
    return os;

}

}
