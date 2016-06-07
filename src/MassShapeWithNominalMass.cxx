#include "MassShapeWithNominalMass.h"

#include "CachedDataValue.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"

namespace yap {

//-------------------------
std::shared_ptr<RealParameter> MassShapeWithNominalMass::mass()
{
    if (!resonance())
        throw exceptions::ResonanceUnset("MassShapeWithNominalMass::mass");
    return resonance()->mass();
}

//-------------------------
void MassShapeWithNominalMass::setParameters(const ParticleTableEntry& entry)
{
    try {
        mass()->setValue(entry.Mass);
    } catch (const exceptions::ResonanceUnset&) { /* ignore */ }
}

//-------------------------
void MassShapeWithNominalMass::setResonance(Resonance* r)
{
    MassShape::setResonance(r);
    addParameter(mass());
}

}




