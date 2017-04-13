#include "Particle.h"

#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
Particle::Particle(const ParticleTableEntry& pde) :
    Particle(pde.name(), pde.quantumNumbers())
{
}

//-------------------------
bool Particle::consistent() const
{
    bool C = true;

    if (ParticleCombinations_.empty()) {
        FLOG(ERROR) << "ParticleCombinations_ is empty";
        C &= false;
    }

    return C;
}

//-------------------------
void Particle::addParticleCombination(const ParticleCombination& pc)
{
    ParticleCombinations_.insert(pc.shared_from_this());
}

//-------------------------
const SpinVector spins(const ParticleVector& v)
{
    SpinVector s;
    s.reserve(v.size());
    std::transform(v.begin(), v.end(), std::back_inserter(s), [](const auto& p){return p->quantumNumbers().twoJ();});
    return s;
}

//-------------------------
const bool decays_to_full_final_state(const Particle& p)
{
    return std::any_of(p.particleCombinations().begin(), p.particleCombinations().end(),
                       [&p](const auto& pc){return pc->indices().size() == p.model()->finalStateParticles().size();});
}

//-------------------------
std::string to_string(const ParticleVector& p, const SpinProjectionVector& two_m)
{
    if (p.size() != two_m.size())
        throw exceptions::Exception("vector size mismatch", "two_string(ParticleVector, SpinProjectionVector)");
    std::string s;
    for (size_t i = 0; i < p.size(); ++i)
        s += ", " + p[i]->name() + " [m = " + spin_to_string(two_m[i]) + "]";
    return s.erase(0, 2);
}

}
