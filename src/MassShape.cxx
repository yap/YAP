#include "MassShape.h"

#include "Exceptions.h"
#include "logging.h"
#include "ParticleCombination.h"
#include "Resonance.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    DataAccessor(&ParticleCombination::equivByOrderlessContent),
    Resonance_(nullptr),
    T_(ComplexCachedDataValue::create(this))
{}

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
void MassShape::addToModel()
{
    DataAccessor::addToModel();
    setDependenciesFromModel();
}

//-------------------------
void MassShape::setResonance(Resonance* r)
{
    if (Resonance_)
        throw exceptions::Exception("MassShape already has owning Resonance", "MassShape::setResonance");

    Resonance_ = r;

    if (Resonance_)
        setDependenciesFromResonance();
}

//-------------------------
void MassShape::replaceResonanceMass(std::shared_ptr<RealParameter> m)
{
    if (!Resonance_)
        throw exceptions::Exception("Resonance unset", "MassShape::replaceResonanceMass");
    Resonance_->setMass(m);
}

//-------------------------
Model* MassShape::model()
{ return (Resonance_) ? Resonance_->model() : nullptr; }

}
