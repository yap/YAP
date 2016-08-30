#include "InvariantMassBinning.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPoint.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include <algorithm>
#include <functional>
#include <vector>

namespace yap {

//-------------------------
InvariantMassBinning::InvariantMassBinning(Model& m, const std::vector<double>& bins) :
    StaticDataAccessor(m, equal_by_orderless_content),
    BinLowEdges_(bins),
    Bin_(RealCachedValue::create(*this))
{
    if (BinLowEdges_.size() < 2)
        throw exceptions::Exception("At least two bin edges must be specified",
                                    "InvariantMassBinning::InvariantMassBinning");

    // use std::less_equal so that if two elements are the same, std::is_sorted
    // will return false
    if (!std::is_sorted(BinLowEdges_.cbegin(), BinLowEdges_.cend(), std::less_equal<double>()))
        throw exceptions::Exception("Elements of the vector are not monotonically increasing",
                                    "InvariantMassBinning::InvariantMassBinning");

    registerWithModel();
}

//-------------------------
void InvariantMassBinning::calculate(DataPoint& d, StatusManager& sm) const
{
    // set all bins to uncalculated
    sm.set(*this, CalculationStatus::uncalculated);

    // loop over particle combinations -> indices
    for (auto& pc_i : symmetrizationIndices()) {
        // check if calculation is necessary
        if (sm.status(*Bin_, pc_i.second ) == CalculationStatus::uncalculated) {
            auto invariant_mass = model()->fourMomenta()->m(d, pc_i.first);

            // get the first lower bound that is greater than the invariant mass
            auto bin = std::upper_bound(BinLowEdges_.cbegin(), BinLowEdges_.cend(), invariant_mass);

            // set the value into the DataPoint
            // NOTE: std::distance() - 1 accounts for the fact that the actual bin value
            // is the one preceding the one that is found above
            Bin_->setValue(std::distance(BinLowEdges_.cbegin(), bin) - 1, d, pc_i.second, sm);
        }
    }
}

//-------------------------
double InvariantMassBinning::bin(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return Bin_->value(d, symmetrizationIndex(pc));
}

}
