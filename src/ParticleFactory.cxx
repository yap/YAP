#include "ParticleFactory.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "MassShape.h"
#include "Resonance.h"

#include <fstream>

namespace yap {

//-------------------------
ParticleTableEntry::ParticleTableEntry(int pdg, std::string name, QuantumNumbers q, double mass, std::vector<double> parameters) :
    QuantumNumbers(q),
    PDG(pdg),
    Name(name),
    Mass(mass),
    MassShapeParameters(parameters)
{
}

//-------------------------
bool ParticleTableEntry::consistent() const
{
    bool C = QuantumNumbers::consistent();

    if (Name.empty()) {
        FLOG(ERROR) << "No name specified.";
        C &= false;
    }

    return C;
}

//-------------------------
std::shared_ptr<FinalStateParticle> ParticleFactory::fsp(int PDG) const
{
    const auto& p = particleTableEntry(PDG);
    DEBUG("make FinalStateParticle " << p.Name << " with quantum numbers " << p);
    return FinalStateParticle::create(p, p.Mass, p.Name);
}

//-------------------------
std::shared_ptr<DecayingParticle> ParticleFactory::decayingParticle(int PDG, double radialSize) const
{
    const auto& p = particleTableEntry(PDG);
    DEBUG("make DecayingParticle " << p.Name << " with quantum numbers " << p);
    return DecayingParticle::create(p, p.Mass, p.Name, radialSize);
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::resonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape) const
{
    const auto& p = particleTableEntry(PDG);
    DEBUG("make Resonance " << p.Name << " with quantum numbers " << p);
    massShape->setParameters(p);
    return Resonance::create(p, p.Mass, p.Name, radialSize, std::move(massShape));
}

//-------------------------
const ParticleTableEntry& ParticleFactory::particleTableEntry(int PDG) const
{
    if (particleTable_.count(PDG) == 0) {
        LOG(FATAL) << "ParticleFactory::particleTableEntry : No particle table entry for PDG " << PDG;
    }
    return particleTable_.at(PDG);
}

//-------------------------
std::pair<ParticleFactory::iterator, bool> ParticleFactory::insert(const value_type& entry)
{
    if (!entry.second.consistent())
        throw exceptions::Exception("entry with PDG code " + std::to_string(entry.first) + " is inconsistent",
                                    "ParticlFactory::insert");

    if (particleTable_.count(entry.first) != 0)
        LOG(WARNING) << "ParticleFactory::insert : PDG code " << entry.first << " already exists. Overwriting entry.";

    return particleTable_.insert(entry);
}

//-------------------------
int ParticleFactory::pdgCode(std::string name) const
{
    auto it = std::find_if(particleTable_.begin(), particleTable_.end(),
    [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.Name == name;});
    if (it == particleTable_.end()) {
        LOG(ERROR) << "ParticleFactory::pdgCode - particle with name \"" << name << "\" not found.";
        return 0;
    }
    return it->first;
}

ParticleFactory read_pdl_file(const std::string& filename)
{
    std::ifstream input(filename.c_str(), std::ios::in);
    return ParticleFactory(PDLIterator(input), PDLIterator::end());
}

}
