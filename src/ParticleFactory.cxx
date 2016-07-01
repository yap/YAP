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
ParticleFactory::ParticleFactory(const std::string pdlFile)
{
    readPDT(pdlFile);
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

//-------------------------
const ParticleTableEntry& ParticleFactory::particleTableEntry(int PDG) const
{
    if (ParticleTable_.count(PDG) == 0)
        throw exceptions::Exception("ParticleFactory::particleTableEntry : No particle table entry for PDG " + std::to_string(PDG),
                                    "ParticleFactory::particleTableEntry");
    return ParticleTable_.at(PDG);
}

//-------------------------
void ParticleFactory::addParticleTableEntry(ParticleTableEntry entry)
{
    if (!entry.consistent())
        throw exceptions::Exception("entry with PDG code " + std::to_string(entry.PDG) + " is inconsistent",
                                    "ParticlFactory::addParticleTableEntry");

    if (ParticleTable_.count(entry.PDG) != 0)
        FLOG(WARNING) << "PDG code " << entry.PDG << " already exists. Overwriting entry.";

    ParticleTable_[entry.PDG] = entry;
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

//-------------------------
void ParticleFactory::readPDT(const std::string fname)
{

    /**
     * This function was taken from EvtGen and modified
     *
     * // Copyright Information: See EvtGen/COPYRIGHT
     * //      Copyright (C) 1998      Caltech, UCSB
     *
     */

    std::ifstream indec;

    indec.open(fname.c_str());

    char cmnd[100];
    char xxxx[100];

    char pname[100];
    int  stdhepid;
    double mass;
    double pwidth;
    double pmaxwidth;
    int    chg3;
    int    spin2;
    double ctau;
    int    lundkc;
    //EvtId i;

    if (!indec)
        throw exceptions::Exception("Could not open \"" + fname + "\"", "ParticleFactory::readPDT");

    do {

        char ch, ch1;

        // ignoring commented lines
        do {
            indec.get(ch);
            if (ch == '\n') {
                indec.get(ch);
            }
            if (ch != '*') {
                indec.putback(ch);
            } else {
                while (indec.get(ch1), ch1 != '\n');
            }
        } while (ch == '*');

        indec >> cmnd;

        if (strcmp(cmnd, "end")) {
            if (!strcmp(cmnd, "add")) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> pname;
                indec >> stdhepid;
                indec >> mass;
                indec >> pwidth;
                indec >> pmaxwidth;
                indec >> chg3;
                indec >> spin2;
                indec >> ctau;
                indec >> lundkc;

                // note: isospin & parity are missing from .pdl format
                addParticleTableEntry(ParticleTableEntry(stdhepid, pname,
                                      QuantumNumbers(spin2, std::round(1. * chg3 / 3)),
                                      mass, {pwidth}));
            }

            // if find a set read information and discard it
            if (!strcmp(cmnd, "set")) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
            }
        }

    } while (strcmp(cmnd, "end"));
}

}
