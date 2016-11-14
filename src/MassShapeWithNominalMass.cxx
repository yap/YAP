#include "MassShapeWithNominalMass.h"

#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
MassShapeWithNominalMass::MassShapeWithNominalMass(double m) :
    MassShape(),
    Mass_(std::make_shared<RealParameter>(m))
{
    addParameter(Mass_);
}

//-------------------------
MassShapeWithNominalMass::MassShapeWithNominalMass(const ParticleTableEntry& pde) :
    MassShapeWithNominalMass(pde.mass())
{
}

}




