#include "MassShapeWithNominalMass.h"

#include "Parameter.h"
#include "ParticleFactory.h"

namespace yap {

//-------------------------
MassShapeWithNominalMass::MassShapeWithNominalMass(double m) :
    MassShape(),
    Mass_(std::make_shared<RealParameter>(m))
{
    addParameter(Mass_);
}

//-------------------------
void MassShapeWithNominalMass::setParameters(const ParticleTableEntry& entry)
{
    if (Mass_->value() < 0)
        *Mass_ = entry.mass();
}

}




