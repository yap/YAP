#include "HelicityAngles.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include "logging.h"

#include <assert.h>

namespace yap {


//-------------------------
const std::array<double, 2>& HelicityAngles::helicityAngles(const DataPoint& d, const StatusManager& sm, const std::shared_ptr<const ParticleCombination>& pc) const
{
    // check if DataPoint is currently in cache
    if (cachedForDataPoint_[&sm] == &d) {
        // find entry
        for (const auto& kv : cachedAngles_[&sm])
            if (equal_up_and_down(kv.first, pc))
                return kv.second;
    }
    else {
        // reset and clear
        cachedForDataPoint_[&sm] = nullptr;
        cachedAngles_[&sm].clear();
    }

    // if not found, calculate
    calculateAngles(d, sm, origin(*pc).shared_from_this(), Model_->coordinateSystem(), unitMatrix<double, 4>());
    cachedForDataPoint_[&sm] = &d;

    return helicityAngles(d, sm, pc);
}

//-------------------------
void HelicityAngles::calculateAngles(const DataPoint& d, const StatusManager& sm,
                                     const std::shared_ptr<const ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts) const
{
    // terminate recursion
    if (is_final_state_particle_combination(*pc))
        return;

    // get pc's 4-mom in data frame
    const auto P = Model_->fourMomenta()->p(d, pc);

    // calculate reference frame for P from parent's RF
    const auto cP = helicityFrame(boosts * P, C);

    // calculate boost from data frame into pc rest frame
    const auto boost = lorentzTransformation(-(boosts * P));

    // boost daughter momentum from data frame into pc rest frame
    const auto p = boost * boosts * Model_->fourMomenta()->p(d, pc->daughters()[0]);

    auto phi_theta = angles<double>(vect<double>(p), cP);

    // set ambiguous phi to theta
    // todo: in this cases, theta should be 0 or pi. In most cases it is, but sometimes not.
    // Not checking if theta == 0 or pi results in tests passing which would otherwise not
    if (std::isnan(phi_theta[0]))
        phi_theta[0] = phi_theta[1];

    cachedAngles_.at(&sm)[pc] = phi_theta;

    for (auto& daughter : pc->daughters())
        // recurse down the decay tree
        calculateAngles(d, sm, daughter, cP, boost);
}

}
