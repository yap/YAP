#include "Particle.h"

#include "logging.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double m, std::string name) :
    AmplitudeComponent(),
    ReportsModel(),
    ReportsParticleCombinations(),
    std::enable_shared_from_this<Particle>(),
    QuantumNumbers_(q),
    Mass_(new RealParameter(m)),
    Name_(name)
{}

//-------------------------
void Particle::setMass(std::shared_ptr<RealParameter> m)
{
    if (!m)
        throw exceptions::Exception("mass is unset", "Particle::setMass");
    Mass_ = m;
}

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
