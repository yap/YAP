#include "MassShape.h"

#include "Constants.h"

#include <algorithm>
#include <set>

namespace yap {

//-------------------------
bool MassShape::areEqual(std::shared_ptr<ParticleCombination> A, std::shared_ptr<ParticleCombination> B) const
{
    // check size of indices vectors
    if (A->indices().size() != B->indices().size())
        return false;

    // check contents of indices vectors
    // (creating a set will sort entries for easy comparison,
    // since order doesn't matter)
    std::set<ParticleIndex> a(A->indices().begin(), A->indices().end());
    std::set<ParticleIndex> b(B->indices().begin(), B->indices().end());
    return std::equal(a.begin(), a.end(), b.begin());
}

}
