#include "Particle.h"
#include "logging.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double mass, std::string name) :
    QuantumNumbers_(q),
    Mass_(mass),
    Name_(name)
{}

//-------------------------
bool Particle::consistent() const
{

    if (Mass_ <= 0.) {
        LOG(ERROR) << "Particle::consistent() - mass not positive.";
        return false;
    }

    return true;
}

}
