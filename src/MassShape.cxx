#include "MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    RecalculableAmplitudeComponent(equal_by_orderless_content),
    Owner_(nullptr),
    T_(ComplexCachedValue::create(*this))
{}

//-------------------------
void MassShape::calculate(DataPartition& D) const
{
    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_si : symmetrizationIndices()) {

        // recalculate & cache, if necessary
        if (D.status(*T(), pc_si.second) == CalculationStatus::uncalculated) {

            calculateT(D, pc_si.first, pc_si.second);
            D.status(*T(), pc_si.second) = CalculationStatus::calculated;
        }

    }
}

//-------------------------
void MassShape::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*T(), CalculationStatus::uncalculated);
}

//-------------------------
const std::complex<double> MassShape::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}

//-------------------------
bool MassShape::consistent() const
{
    bool C = DataAccessor::consistent();

    // check if owner is set
    if (!Owner_) {
        FLOG(ERROR) << "Owner isn't set";
        C &= false;
    }

    return C;
}

//-------------------------
void MassShape::setOwner(DecayingParticle* dp)
{
    if (Owner_)
        throw exceptions::Exception("MassShape already has owner", "MassShape::setOwner");

    Owner_ = dp;
}

//-------------------------
const Model* MassShape::model() const
{
    return Owner_->model();
}

}
