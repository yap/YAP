#include "BreitWigner.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double m, double w) :
    MassShapeWithNominalMass(m),
    T_(ComplexCachedValue::create(*this)),
    Width_(std::make_shared<PositiveRealParameter>(w))
{
    addParameter(Width_);
}

//-------------------------
BreitWigner::BreitWigner(const ParticleTableEntry& pde) :
    BreitWigner(pde.mass(), get_nth_element(pde, 0, "BreitWigner::BreitWigner"))
{
}

//-------------------------
void BreitWigner::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*T_, CalculationStatus::uncalculated);
}

//-------------------------
const std::complex<double> BreitWigner::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void BreitWigner::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T_, si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:
    
    // mass * width
    auto mw = mass()->value() * Width_->value();

    // mass^2 - i * mass * width
    auto m2_imw = pow(mass()->value(), 2) - 1_i * mass()->value() * Width_->value();

    /////////////////////////
    
    // T := mass * width / (mass^2 - s - i * mass * width)
    for (auto& d : D)
        T_->setValue(mw / (m2_imw - model()->fourMomenta()->m2(d, pc)), d, si, D);

    D.status(*T_, si) = CalculationStatus::calculated;
}

}




