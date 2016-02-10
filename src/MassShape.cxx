#include "MassShape.h"

#include "Exceptions.h"
#include "DecayingParticle.h"
#include "Resonance.h"

namespace yap {

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

    if (Resonance_)
        borrowParametersFromResonance();
}

//-------------------------
Model* MassShape::model()
{ return (Resonance_) ? Resonance_->model() : nullptr; }

}
