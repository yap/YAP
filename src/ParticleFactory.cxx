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
ParticleFactory::ParticleFactory(const std::string pdlFile) :
    InitialStateParticle_(nullptr)
{
    readPDT(pdlFile);
}

//-------------------------
std::shared_ptr<FinalStateParticle> ParticleFactory::createFinalStateParticle(int PDG, std::vector<ParticleIndex> indices)
{
    const PdlParticleProperties& p = particleProperties(PDG);
    return std::make_shared<FinalStateParticle>(createQuantumNumbers(PDG), p.Mass_, p.Name_, PDG, indices);
}

//-------------------------
std::shared_ptr<InitialStateParticle> ParticleFactory::createInitialStateParticle(int PDG, double radialSize)
{
    const PdlParticleProperties& p = particleProperties(PDG);
    std::shared_ptr<InitialStateParticle> isp = std::make_shared<InitialStateParticle>(createQuantumNumbers(PDG), p.Mass_, p.Name_, radialSize);
    InitialStateParticle_ = isp.get();

    return isp;
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::createResonance(int PDG, double radialSize, MassShape* massShape)
{
    const PdlParticleProperties& p = particleProperties(PDG);
    return std::make_shared<Resonance>(createQuantumNumbers(PDG), p.Mass_, p.Name_, radialSize, massShape);
}

//-------------------------
std::shared_ptr<Resonance> ParticleFactory::createResonanceBreitWigner(int PDG, double radialSize)
{
    const PdlParticleProperties& p = particleProperties(PDG);
    return std::make_shared<Resonance>( createQuantumNumbers(PDG),
                                        p.Mass_, p.Name_, radialSize,
                                        new BreitWigner(initialStateParticle(), p.Mass_, p.Width_) );
}

//-------------------------
QuantumNumbers ParticleFactory::createQuantumNumbers(int PDG)
{
    const PdlParticleProperties& p = particleProperties(PDG);

    unsigned char twoJ = p.TwoJ_;
    char Q = std::round(1. / 3. * p.ThreeCharge_);

    // \todo
    unsigned char twoI = 0;
    unsigned char P = 0;

    return QuantumNumbers(twoI, twoJ, P, Q);
}

//-------------------------
const PdlParticleProperties& ParticleFactory::particleProperties(int PDG) const
{
    if (particleProperties_.count(PDG) == 0) {
        LOG(ERROR) << "No PdlParticleProperties for PDG " << PDG;
    }
    return particleProperties_.at(PDG);
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

        do {

            indec.get(ch);
            if (ch == '\n') indec.get(ch);
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


                PdlParticleProperties tmp;

                tmp.PDGCode_ = stdhepid;
                tmp.Name_ = std::string(pname);
                tmp.Mass_ = mass;
                tmp.Width_ = pwidth;
                tmp.ThreeCharge_ = chg3;
                tmp.TwoJ_ = spin2;

                particleProperties_[stdhepid] = tmp;

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
