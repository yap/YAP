#include "MassShape.h"

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
    Resonance_ = r;

    if (Resonance_)
        borrowParametersFromResonance();
}

//-------------------------
InitialStateParticle* MassShape::initialStateParticle()
{ return (Resonance_) ? Resonance_->initialStateParticle() : nullptr; }

}
