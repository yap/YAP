#include "BreitWigner.h"

#include "CachedValue.h"
#include "Constants.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double w) :
    MassShapeWithNominalMass(),
    Width_(std::make_shared<RealParameter>(w))
{
    addParameter(Width_);
}

//-------------------------
void BreitWigner::setParameters(const ParticleTableEntry& entry)
{
    MassShapeWithNominalMass::setParameters(entry);

    if (entry.MassShapeParameters.empty())
        throw exceptions::Exception("entry.MassShapeParameter is empty", "BreitWigner::setParameters");

    if (width()->value() < 0)
        width()->setValue(entry.MassShapeParameters[0]);
}

//-------------------------
void BreitWigner::calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const
{
    // common factor := M^2 - i * M * Gamma
    auto M2_iMG = pow(mass()->value(), 2) - Complex_i * mass()->value() * width()->value();

    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D)
        T()->setValue(1. / (M2_iMG - model()->fourMomenta()->m2(d, pc)), d, si, D);
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool C = MassShapeWithNominalMass::consistent();

    if (width()->value() <= 0) {
        FLOG(ERROR) << "width <= 0";
        C &= false;
    }

    return C;
}

}




