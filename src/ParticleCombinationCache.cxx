#include "ParticleCombinationCache.h"

#include "container_utils.h"
#include "Exceptions.h"
#include "logging.h"

#include <algorithm>
#include <iomanip>

namespace yap {

//-------------------------
/// \todo Shouldn't this management be done in the creation process?
ParticleCombinationCache::ParticleCombinationCache(std::vector<shared_ptr_type> V)
    : WeakPtrCache(V)
{
    // // set parents from isp down
    // for (auto& pc : V)
    //     setLineage(pc);

    // check consistency
    if (!consistent())
        throw exceptions::Exception("Inconsistent ParticleCombinationCache", "ParticleCombinationCache::ParticleCombinationCache");

    // remove expired cache entries
    removeExpired();
}

//-------------------------
ParticleCombinationCache::shared_ptr_type ParticleCombinationCache::create_copy(const ParticleCombination& other) const
{
    auto pc = shared_ptr_type(new ParticleCombination());

    // if copying a final-state particle, just copy indices
    if (other.isFinalStateParticle())
        pc->Indices_ = other.Indices_;

    // else copy daughters (recursively) and add them
    else
        for (auto& d : other.Daughters_)
            pc->addDaughter(create_copy(*d));

    return pc;
}

//-------------------------
ParticleCombinationCache::shared_ptr_type ParticleCombinationCache::create_composite(const ParticleCombinationVector& D) const
{
    auto pc = shared_ptr_type(new ParticleCombination());

    for (auto& d : D)
        pc->addDaughter(create_copy(*d));

    return pc;
}

//-------------------------
// void ParticleCombinationCache::setLineage(shared_ptr_type pc)
// {
//     /// \todo why do we need to do this here?
//     for (auto& D : pc->Daughters_) {
//         // copy D
//         auto d = std::make_shared<type>(*D);
//         // set copy's parent to pc
//         std::const_pointer_cast<ParticleCombination>(d)->Parent_ = pc;
//         // call recursively
//         setLineage(pc);
//         // replace daughter with copy (from cache)
//         std::const_pointer_cast<ParticleCombination>(D) = std::const_pointer_cast<ParticleCombination>(operator[](d));
//         // std::const_pointer_cast<ParticleCombination>(D).swap(std::const_pointer_cast<ParticleCombination>(operator[](d)));
//     }
// }

//-------------------------
ParticleCombinationCache::weak_ptr_type ParticleCombinationCache::findByUnorderedContent(const std::vector<ParticleIndex>& I) const
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
    for (auto& d : pc->daughters())
        addToCache(d);
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
std::ostream& print_daughters(std::ostream& os, const ParticleCombinationCache& C, std::shared_ptr<ParticleCombination> pc, unsigned ndig, std::string prefix, std::set<unsigned>& used)
{
    for (auto& d : pc->daughters()) {
        unsigned i = 0;
        for (auto& w : C) {
            if (!w.expired() and w.lock() == d) {
                os << std::setfill('0') << std::setw(ndig) << i << " : " << prefix << *w.lock();
                if (d->parent() != pc)
                    os << " (parent unset)";
                os << std::endl;
                used.insert(i);
                break;
            }
            ++i;
        }
        print_daughters(os, C, d, ndig, prefix + "    ", used);
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


    // print final state particles
    // unsigned i = 0;
    // for (auto& w : *this) {
    //     if (!w.expired() and w.lock()->isFinalStateParticle())
    //         os << std::setfill('0') << std::setw(ndig) << i << " : (fsp) " << *w.lock() << std::endl;
    //     ++i;
    // }

    std::set<unsigned> used;

    // print others as decay trees starting from ISP's
    unsigned i = 0;
    for (auto& w : *this) {
        if (!w.expired() and w.lock()->indices().size() == n_fsp) {
            // print isp
            os << std::setfill('0') << std::setw(ndig) << i << " : " << *w.lock() << std::endl;
            used.insert(i);

            // print daughters
            print_daughters(os, *this, w.lock(), ndig, "    ", used);
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

    // for (auto wpc : C) {
    //     auto pc = wpc.lock();
    //     s += to_string(*pc);
    //     auto pt = pc->parent();
    //     if (pt) {
    //         s += " in decay chain ";
    //         while (pt->parent())
    //             pt = pt->parent();
    //         s += to_string(*pt);
    //     }
    //     s += "\n";
    // }
    // s.erase(s.size() - 1, 1);
    // return s;
}

}
