#include "Particle.h"

#include "logging.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double m, std::string name) :
    AmplitudeComponent(),
    BelongsToInitialStateParticle(),
    QuantumNumbers_(q),
    Mass_(new RealParameter(m)),
    Name_(name)
{}

//-------------------------
bool Particle::consistent() const
{
    if (Mass_->value() <= 0.) {
        LOG(ERROR) << "Particle::consistent() - mass not positive.";
        return false;
    }

    return true;
}

}
