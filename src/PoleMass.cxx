#include "PoleMass.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
PoleMass::PoleMass(std::complex<double> mass) :
    MassShape(),
    T_(ComplexCachedValue::create(*this)),
    Mass_(std::make_shared<ComplexParameter>(mass))
{
    addParameter(Mass_);
}

//-------------------------
PoleMass::PoleMass(const ParticleTableEntry& pde) :
    PoleMass(std::complex<double>(pde.mass(), get_nth_element(pde, 0, "PoleMass::PoleMass") / 2))
{
}

//-------------------------
void PoleMass::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*T_, CalculationStatus::uncalculated);
}

//-------------------------
const std::complex<double> PoleMass::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void PoleMass::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T_, si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    // 2 * re(mass) * im(mass)
    auto two_re_im_m = 2. * real(Mass_->value()) * imag(Mass_->value());

    // mass^2
    auto m2 = pow(Mass_->value(), 2);
    
    /////////////////////////
    
    // T := 2 re(mass) im(mass) / (mass^2 - s)
    for (auto& d : D)
        T_->setValue(two_re_im_m / (m2 - model()->fourMomenta()->m2(d, pc)), d, si, D);

    D.status(*T_, si) = CalculationStatus::calculated;
}

//-------------------------
bool PoleMass::consistent() const
{
    bool C = MassShape::consistent();

    if (real(Mass_->value()) < 0) {
        FLOG(ERROR) << "real(mass) < 0";
        C &= false;
    }
    if (imag(Mass_->value()) >= 0) {
        FLOG(ERROR) << "imag(mass) >= 0";
        C &= false;
    }

    return C;
}

}




