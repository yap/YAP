#include "Particle.h"

#include "logging.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double mass, std::string name) :
    QuantumNumbers_(q),
    Mass_(new RealParameter(mass)),
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
