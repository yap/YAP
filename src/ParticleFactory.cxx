#include "ParticleFactory.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "MassShape.h"
#include "Resonance.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <utility>

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
    return FinalStateParticle::create(p, p.Mass, p.Name);
}

//-------------------------
std::shared_ptr<DecayingParticle> ParticleFactory::decayingParticle(int PDG, double radialSize) const
{
    const auto& p = particleTableEntry(PDG);
    return DecayingParticle::create(p, p.Mass, p.Name, radialSize);
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::resonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape) const
{
    const auto& p = particleTableEntry(PDG);
    massShape->setParameters(p);
    return Resonance::create(p, p.Mass, p.Name, radialSize, std::move(massShape));
}

////-------------------------
//const ParticleFactory& ParticleFactory::operator+=(const ParticleFactory& rhs)
//{
//	std::transform(rhs.ParticleTable_.begin(), rhs.ParticleTable_.end(),
//			inserter(*this), std::bind(std::get<1>, std::placeholders::_1));
//    return *this;
//}

//-------------------------
const ParticleTableEntry& ParticleFactory::particleTableEntry(int PDG) const
{
    if (ParticleTable_.count(PDG) == 0)
        throw exceptions::Exception("ParticleFactory::particleTableEntry : No particle table entry for PDG " + std::to_string(PDG),
                                    "ParticleFactory::particleTableEntry");
    return ParticleTable_.at(PDG);
}

//-------------------------
std::pair<ParticleTableMap::iterator, bool> ParticleFactory::insert(const ParticleTableEntry& entry)
{
    if (!entry.consistent())
        throw exceptions::Exception("entry with PDG code " + std::to_string(entry.PDG) + " is inconsistent",
                                    "ParticlFactory::insert");

    if (ParticleTable_.count(entry.PDG) != 0)
        LOG(WARNING) << "ParticleFactory::insert : PDG code " << entry.PDG << " already exists. Overwriting entry.";

    return ParticleTable_.insert(ParticleTableMap::value_type(entry.PDG, entry));
}

//-------------------------
ParticleTableMap::iterator ParticleFactory::insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry)
{
    if (!entry.consistent())
        throw exceptions::Exception("entry with PDG code " + std::to_string(entry.PDG) + " is inconsistent",
                                    "ParticlFactory::insert");

    if (ParticleTable_.count(entry.PDG) != 0)
        LOG(WARNING) << "ParticleFactory::insert : PDG code " << entry.PDG << " already exists. Overwriting entry.";

    return ParticleTable_.insert(hint, ParticleTableMap::value_type(entry.PDG, entry));
}

//-------------------------
int ParticleFactory::pdgCode(std::string name) const
{
    auto it = std::find_if(ParticleTable_.begin(), ParticleTable_.end(),
    [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.Name == name;});
    if (it == ParticleTable_.end())
        throw exceptions::Exception("particle with name \"" + name + "\" not found",
                                    "ParticleFactory::pdgCode");
    return it->first;
}

}
