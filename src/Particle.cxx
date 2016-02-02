#include "Particle.h"

#include "logging.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double m, std::string name) :
    AmplitudeComponent(),
    ReportsInitialStateParticle(),
    ReportsParticleCombinations(),
    QuantumNumbers_(q),
    Mass_(new RealParameter(m)),
    Name_(name)
{}

//-------------------------
bool Particle::consistent() const
{
    if (Mass_->value() < 0.) {
        FLOG(ERROR) << "mass is negative";
        return false;
    }

    return true;
}

}
