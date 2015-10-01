#include "ParticleFactory.h"

#include "BreitWigner.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <fstream>
#include <iostream>

namespace yap {

//-------------------------
ParticleFactory::ParticleTableEntry::ParticleTableEntry(int pdg, std::string name, QuantumNumbers q, double mass, double width) :
    QuantumNumbers(q),
    PDG_(pdg),
    Name_(name),
    Mass_(mass),
    Width_(width)
{
}

//-------------------------
bool ParticleFactory::ParticleTableEntry::consistent() const
{
    bool result = QuantumNumbers::consistent();

    if (Name_.empty()) {
        LOG(ERROR) << "ParticleTableEntry::consistent : No name specified.";
        result = false;
    }

    return result;
}

//-------------------------
ParticleFactory::ParticleFactory(const std::string pdlFile) :
    InitialStateParticle_(nullptr)
{
    readPDT(pdlFile);
}

//-------------------------
std::shared_ptr<FinalStateParticle> ParticleFactory::createFinalStateParticle(int PDG, std::vector<ParticleIndex> indices)
{
    const ParticleTableEntry& p = particleTableEntry(PDG);
    DEBUG("make FinalStateParticle " << p.Name_ << " with quantum numbers " << p);
    return std::make_shared<FinalStateParticle>(p, p.Mass_, p.Name_, indices);
}

//-------------------------
std::shared_ptr<InitialStateParticle> ParticleFactory::createInitialStateParticle(int PDG, double radialSize)
{
    const ParticleTableEntry& p = particleTableEntry(PDG);

    if (p.twoJ() != 0)
        LOG(ERROR) << "InitialStateParticle has spin != 0. ";

    DEBUG("make InitialStateParticle " << p.Name_ << " with quantum numbers " << p);
    std::shared_ptr<InitialStateParticle> isp = std::make_shared<InitialStateParticle>(p, p.Mass_, p.Name_, radialSize);
    InitialStateParticle_ = isp.get();

    return isp;
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::createResonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape)
{
    const ParticleTableEntry& p = particleTableEntry(PDG);
    DEBUG("make Resonance " << p.Name_ << " with quantum numbers " << p);
    return std::make_shared<Resonance>(p, p.Mass_, p.Name_, radialSize, massShape);
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::createResonanceBreitWigner(int PDG, double radialSize)
{
    const ParticleTableEntry& p = particleTableEntry(PDG);
    std::shared_ptr<MassShape> massShape = std::make_shared<BreitWigner>(initialStateParticle(), p.Mass_, p.Width_);
    return createResonance(PDG, radialSize, massShape);
}

//-------------------------
const ParticleFactory::ParticleTableEntry& ParticleFactory::particleTableEntry(int PDG) const
{
    if (particleTable_.count(PDG) == 0) {
        LOG(FATAL) << "ParticleFactory::particleTableEntry : No particle table entry for PDG " << PDG;
    }
    return particleTable_.at(PDG);
}

//-------------------------
bool ParticleFactory::addParticleTableEntry(ParticleTableEntry entry)
{
    if (!entry.consistent()) {
        LOG(ERROR) << "ParticleFactory::addParticleTableEntry : entry with PDG code " << entry.PDG_ << " inconsistent";
        return false;
    }

    if (particleTable_.count(entry.PDG_) != 0) {
        LOG(WARNING) << "ParticleFactory::addParticleTableEntry : PDG code " << entry.PDG_ << " already exists. Overwriting entry.";
    }

    particleTable_[entry.PDG_] = entry;
    // } else
    //     particleTable_.insert(std::make_pair<int, ParticleTableEntry>(entry.PDG_, entry));

    return true;
}

//-------------------------
InitialStateParticle* ParticleFactory::initialStateParticle()
{
    if (InitialStateParticle_ == nullptr)
        LOG(ERROR) << "ParticleFactory: no initialStateParticle. createInitialStateParticle first before creating other resonances.";

    return InitialStateParticle_;
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

    if (!indec) {
        LOG(ERROR) << "Could not open:" << fname.c_str() << "EvtPDL";
        return;
    }

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
                                      QuantumNumbers(0, spin2, 0, std::round(1. * chg3 / 3)),
                                      mass, pwidth));
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
