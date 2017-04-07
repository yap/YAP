#include "HelicityAngles.h"

#include "FourMomenta.h"
#include "FourVector.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

namespace yap {


//-------------------------
const spherical_angles<double>& HelicityAngles::operator()(const DataPoint& d, const StatusManager& sm, const std::shared_ptr<const ParticleCombination>& pc) const
{
    // check if StatusManager is in cache
    auto it = CachedAngles_.find(&sm);
    if (CachedAngles_.find(&sm) == CachedAngles_.end()) {
        addToCache(sm);
        it = CachedAngles_.find(&sm);
    }

    // get entries for status manager
    angles_cache& cachedAngles = it->second;
    auto& cachedForDataPoint = CachedForDataPoint_.at(&sm);


    // check if DataPoint is currently in cache
    if (cachedForDataPoint == &d) {
        // find and return entry
        for (const auto& pc_angles : cachedAngles)
            if (equal_up_and_down(pc_angles.first, pc))
                return pc_angles.second;
    }
    else {
        // reset and clear
        cachedForDataPoint = nullptr;
        cachedAngles.clear();
    }

    // if not found, calculate
    calculateAngles(d, sm, origin(*pc).shared_from_this(), Model_->coordinateSystem(), unitMatrix<double, 4>());
    cachedForDataPoint = &d;

    // find and return entry
    for (const auto& pc_angles : cachedAngles)
        if (equal_up_and_down(pc_angles.first, pc))
            return pc_angles.second;

    throw exceptions::Exception("Cannot find entry for DataPoint and ParticleCombination even though it should have been calculated. Something went wrong!", "HelicityAngles::operator()");
}

//-------------------------
void HelicityAngles::addToCache(const StatusManager& sm) const
{
    std::lock_guard<std::mutex> guard(CacheMutex_);
    CachedForDataPoint_[&sm] = nullptr;
    CachedAngles_[&sm].clear();
}

//-------------------------
void HelicityAngles::clearCache(const StatusManager& sm) const 
{
    std::lock_guard<std::mutex> guard(CacheMutex_);
    if (CachedAngles_.find(&sm) != CachedAngles_.end()) {
        CachedAngles_.at(&sm).clear();
        CachedForDataPoint_.at(&sm) = nullptr;
    }
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

    auto hel_angles = angles<double>(vect<double>(p), cP);

    // set ambiguous phi to theta
    // todo: in this cases, theta should be 0 or pi. In most cases it is, but sometimes not.
    // Not checking if theta == 0 or pi results in tests passing which would otherwise not
    if (std::isnan(hel_angles.phi))
        hel_angles.phi = hel_angles.theta;

    CachedAngles_.at(&sm)[pc] = hel_angles;

    for (auto& daughter : pc->daughters())
        // recurse down the decay tree
        calculateAngles(d, sm, daughter, cP, boost);
}

}
