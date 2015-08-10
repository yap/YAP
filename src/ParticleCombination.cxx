#include "ParticleCombination.h"

#include <algorithm>

namespace yap {

//-------------------------
    bool ParticleCombination::consistent() const
    {
        // create unique_copy of this object
        ParticleCombination U;
        std::unique_copy(begin(), end(), U.begin());
        // check unique_copy-object's size == this object's size
        return U.size() == size();
    }
    
    //-------------------------
    bool operator==(const ParticleCombination& A, const ParticleCombination& B)
    {
        if (A.size() != B.size())
            return false;
        return std::equal(A.begin(), A.end(), B.begin());
    }

}
