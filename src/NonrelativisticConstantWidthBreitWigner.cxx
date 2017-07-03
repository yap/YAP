#include "NonrelativisticConstantWidthBreitWigner.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void NonrelativisticConstantWidthBreitWigner::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:
    
    // mass - i  width
    auto m_iw = mass()->value() - 1_i * width()->value() / 2.;

    /////////////////////////
    
    // T := width / (mass - sqrt(s) - i * width)
    for (auto& d : D)
        T()->setValue(width()->value() / 2. / (m_iw - model()->fourMomenta()->m(d, pc)), d, si, D);

    D.status(*T(), si) = CalculationStatus::calculated;
}

}




