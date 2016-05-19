#include "ReportsParticleCombinations.h"

#include <algorithm>

namespace yap {

//-------------------------
bool ReportsParticleCombinations::hasParticleCombination(const std::shared_ptr<ParticleCombination>& PC,
        const ParticleCombination::Equiv& equiv) const
{
    const auto& pcv = particleCombinations();
    return std::find_if(pcv.begin(), pcv.end(),
                        [&](const ParticleCombinationVector::value_type & pc)
    { return equiv(PC, pc); }) != pcv.end();
}

}
