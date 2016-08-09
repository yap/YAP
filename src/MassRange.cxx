#include "MassRange.h"

#include "container_utils.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "MassAxes.h"
#include "Parameter.h"
#include "ParticleCombination.h"

#include <algorithm>
#include <functional>

namespace yap {

//-------------------------
const FinalStateParticleVector::value_type& fsp_with_index(const FinalStateParticleVector& FSPs, unsigned i)
{
    for (const auto& fsp : FSPs)
        for (const auto& pc : fsp->particleCombinations())
            if (pc->indices()[0] == i)
                return fsp;
    throw exceptions::Exception("fsp with index " + std::to_string(i) + " not found", "fsp_with_index");
}

//-------------------------
const MassRange mass_range(const std::shared_ptr<ParticleCombination>& pc, std::shared_ptr<DecayingParticle> ISP, const FinalStateParticleVector& FSPs)
{
    // get list of PCs of ISP that decay to FSPs
    ParticleCombinationVector isp_pcs;
    isp_pcs.reserve(ISP->particleCombinations().size());
    std::copy_if(ISP->particleCombinations().begin(), ISP->particleCombinations().end(),
                 std::back_inserter(isp_pcs), std::bind(valid_final_state, std::placeholders::_1, FSPs));
    if (isp_pcs.empty())
        throw exceptions::Exception("FSPs not a final state of ISP", "mass_range");

    // pc must be a subset of one of these PCs
    auto it = std::find_if(isp_pcs.begin(), isp_pcs.end(), [&](const std::shared_ptr<ParticleCombination>& isp) {return contains(isp->indices(), pc->indices());});
    if (it == isp_pcs.end())
        throw exceptions::Exception("pc does not match ISP and FPSs", "mass_range");

    std::array<double, 2> m = {0, ISP->mass()->value()};

    for (const auto i : (*it)->indices()) {
        // get fsp
        auto fsp = fsp_with_index(FSPs, i);
        if (std::find(pc->indices().begin(), pc->indices().end(), i) != pc->indices().end())
            // add mass to low end
            m[0] += fsp->mass()->value();
        else
            // subtract mass from high end
            m[1] -= fsp->mass()->value();
    }
    return m;
}

//-------------------------
const std::vector<MassRange> mass_range(const MassAxes& A, std::shared_ptr<DecayingParticle> ISP, const FinalStateParticleVector& FSPs)
{
    std::vector<MassRange> R;
    R.reserve(A.size());
    std::transform(A.begin(), A.end(), std::back_inserter(R), [&](const MassAxes::value_type & a) {return mass_range(a, ISP, FSPs);});
    return R;
}

//-------------------------
const MassRange squared(MassRange mr)
{
    std::transform(mr.begin(), mr.end(), mr.begin(), std::bind(pow, std::placeholders::_1, 2));
    return mr;
}

//-------------------------
const std::vector<MassRange> squared(std::vector<MassRange> mr)
{
    std::transform(mr.begin(), mr.end(), mr.begin(), [](const MassRange & r) {return squared(r);});
    return mr;
}

}
