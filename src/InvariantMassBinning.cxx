#include "InvariantMassBinning.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include <algorithm>
#include <functional>

namespace yap {

//-------------------------
InvariantMassBinning::InvariantMassBinning(Model& m, const std::vector<double>& bins) :
    StaticDataAccessor(m, equal_by_orderless_content),
    Bins_(bins),
    BinNumber_(RealCachedValue::create(*this))
{
    if (Bins_.size() < 2)
        throw exceptions::Exception("At least two bins must be specified",
                                    "InvariantMassBinning::InvariantMassBinning");

    // use std::less_equal so that if two elements are the same, std::is_sorted
    // will return false
    if (!std::is_sorted(Bins_.cbegin(), Bins_.cend(), std::less_equal<decltype(*Bins_.cbegin())>()))
        throw exceptions::Exception("Elements of the vector are not monotonically increasing",
                                    "InvariantMassBinning::InvariantMassBinning");
}

//-------------------------
void InvariantMassBinning::calculate(DataPoint& d, StatusManager& sm) const
{
    // set all bins to uncalculated
    sm.set(*BinNumber_, CalculationStatus::uncalculated);

    // loop over particle combinations -> indices
    for (auto& pc_i : symmetrizationIndices()) {
        // check if calculation is necessary
        if (sm.status(*BinNumber_, pc_i.second ) == CalculationStatus::uncalculated) {
            auto invariant_mass = model()->fourMomenta()->m(d, pc_i.first);

            // get the first lower bound that is greater than the invariant mass
            auto bin = std::upper_bound(Bins_.cbegin(), Bins_.cend(), invariant_mass);

            // check if the value above the upper limit
            if (bin == Bins_.cend())
                throw exceptions::Exception("Mass is above the upper bin",
                                            "InvariantMassBinning::calculate");
            // check if the value below the lower limit
            if (bin == Bins_.cbegin())
                throw exceptions::Exception("Mass is below the lower bin",
                                            "InvariantMassBinning::calculate");

            // if the previous checks passed, set bin to the actual value
            // of the containing bin (i.e. bin is the last value in the vector
            // that is not greater than invariant_mass);
            --bin;

            // set the value into the DataPoint
            BinNumber_->setValue(std::distance(Bins_.cbegin(), bin), d, pc_i.second, sm);
        }
    }
}

}
