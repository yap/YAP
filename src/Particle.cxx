#include "Particle.h"

#include "logging.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double m, std::string name) :
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
    bool C = true;

    if (Mass_->value() < 0.) {
        FLOG(ERROR) << "mass is negative";
        C &= false;
    }

    if (ParticleCombinations_.empty()) {
        FLOG(ERROR) << "ParticleCombinations_ is empty";
        C &= false;
    }

    return C;
}

//-------------------------
void Particle::addParticleCombination(const std::shared_ptr<ParticleCombination>& pc)
{
    ParticleCombinations_.push_back(pc);
}

//-------------------------
void Particle::pruneParticleCombinations()
{
    prune_particle_combinations(ParticleCombinations_);
}

//-------------------------
const SpinVector spins(const ParticleVector& v)
{
    SpinVector s;
    s.reserve(v.size());
    std::transform(v.begin(), v.end(), std::back_inserter(s),
    [](const ParticleVector::value_type & p) {return p->quantumNumbers().twoJ();});
    return s;
}

//-------------------------
std::string to_string(const Particle& p)
{
    return p.name() + "(" + to_string(p.quantumNumbers()) + "), mass = " + std::to_string(p.mass()->value());
}

//-------------------------
const bool orderless_equal(const ParticleVector& A, const ParticleVector& B)
{
    if (A == B)
        return true;
    if (A.size() != B.size())
        return false;
    // create set
    std::set<ParticleVector::value_type> S(A.begin(), A.end());
    return std::all_of(S.begin(), S.end(), [&](const ParticleVector::value_type& p)
                       {return std::count(A.begin(), A.end(), p) == std::count(B.begin(), B.end(), p);});
}

}
