#include "ParticleFactory.h"

#include "DecayingParticle.h"
#include "Exceptions.h"
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
    PDG_(pdg),
    Name_(name),
    QuantumNumbers_(q),
    Mass_(mass),
    MassShapeParameters_(parameters)
{
    if (Name_.empty())
        throw exceptions::Exception("Name is empty", "ParticleTableEntry::ParticleTableEntry");
    if (Mass_ < 0)
        throw exceptions::Exception("Mass is negative", "ParticleTableEntry::ParticleTableEntry");
}

//-------------------------
std::shared_ptr<FinalStateParticle> ParticleFactory::fsp(int PDG) const
{
    const auto& p = (*this)[PDG];
    return FinalStateParticle::create(p.name(), p.quantumNumbers(), p.mass());
}

//-------------------------
std::shared_ptr<DecayingParticle> ParticleFactory::decayingParticle(int PDG, double radialSize) const
{
    const auto& p = (*this)[PDG];
    return DecayingParticle::create(p.name(), p.quantumNumbers(), radialSize);
}

//-------------------------
    std::shared_ptr<Resonance> ParticleFactory::resonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape) const
{
    const auto& p = (*this)[PDG];
    massShape->setParameters(p);
    return Resonance::create(p.name(), p.quantumNumbers(), radialSize, std::move(massShape));
}

//-------------------------
ParticleFactory& ParticleFactory::operator+=(const ParticleFactory& rhs)
{
    std::transform(rhs.ParticleTable_.begin(), rhs.ParticleTable_.end(), inserter(*this),
                   std::bind(&ParticleTableMap::value_type::second, std::placeholders::_1));
    return *this;
}

//-------------------------
const ParticleTableEntry& ParticleFactory::operator[](int PDG) const
{
    if (ParticleTable_.count(PDG) == 0)
        throw exceptions::Exception("No particle table entry for PDG " + std::to_string(PDG), "ParticleFactory::operator[]");
    return ParticleTable_.at(PDG);
}

//-------------------------
std::pair<ParticleTableMap::iterator, bool> ParticleFactory::insert(const ParticleTableEntry& entry)
{
    auto it_b = ParticleTable_.insert(ParticleTableMap::value_type(entry.pdg(), entry));

    // if insertion failed because key value entry.PDG was already contained
    if (!it_b.second and it_b.first != ParticleTable_.end()) {
        LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
        it_b.first->second = entry;
        it_b.second = true;
    }

    return it_b;
}

//-------------------------
ParticleTableMap::iterator ParticleFactory::insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry)
{
    if (ParticleTable_.count(entry.pdg()) == 0)
        return ParticleTable_.insert(hint, ParticleTableMap::value_type(entry.pdg(), entry));

    LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
    auto it = ParticleTable_.find(entry.pdg());
    it->second = entry;
    return it;
}

//-------------------------
int ParticleFactory::pdgCode(std::string name) const
{
    auto it = std::find_if(ParticleTable_.begin(), ParticleTable_.end(),
    [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.name() == name;});
    if (it == ParticleTable_.end())
        throw exceptions::Exception("particle with name \"" + name + "\" not found", "ParticleFactory::pdgCode");
    return it->first;
}

}
