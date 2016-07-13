#include "MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "Resonance.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    RecalculableDataAccessor(equal_by_orderless_content),
    Resonance_(nullptr),
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
    if (variable_status(*this) == VariableStatus::changed)
        D.set(*T(), CalculationStatus::uncalculated);
}

//-------------------------
std::complex<double> MassShape::value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}

//-------------------------
bool MassShape::consistent() const
{
    bool C = DataAccessor::consistent();

    // check if owning resonance is set
    if (!Resonance_) {
        FLOG(ERROR) << "Owning resonance isn't set";
        C &= false;
    }

    return C;
}

//-------------------------
void MassShape::setResonance(Resonance* r)
{
    if (Resonance_)
        throw exceptions::Exception("MassShape already has owning Resonance", "MassShape::setResonance");

    Resonance_ = r;
}

//-------------------------
void MassShape::replaceResonanceMass(std::shared_ptr<RealParameter> m)
{
    if (!Resonance_)
        throw exceptions::Exception("Resonance unset", "MassShape::replaceResonanceMass");
    Resonance_->setMass(m);
}

//-------------------------
const Model* MassShape::model() const
{
    return Resonance_->model();
}

}
