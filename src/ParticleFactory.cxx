#include "ParticleFactory.h"

#include "BreitWigner.h"

#include <TDatabasePDG.h>

#include "logging.h"

namespace yap {

//-------------------------
FinalStateParticle* ParticleFactory::createFinalStateParticle(int PDG)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new FinalStateParticle(createQuantumNumbers(PDG), p->Mass(), std::string(p->GetName()), PDG);
}

//-------------------------
InitialStateParticle* ParticleFactory::createInitialStateParticle(int PDG, double radialSize)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new InitialStateParticle(createQuantumNumbers(PDG), p->Mass(), p->GetName(), radialSize);
}

//-------------------------
Resonance* ParticleFactory::createResonance(int PDG, double radialSize, MassShape* massShape)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new Resonance(createQuantumNumbers(PDG), p->Mass(), p->GetName(), radialSize, massShape);
}

//-------------------------
Resonance* ParticleFactory::createResonanceBreitWigner(int PDG, double radialSize)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new Resonance(createQuantumNumbers(PDG), p->Mass(), std::string(p->GetName()),
                         radialSize, new BreitWigner(p->Mass(), p->Width()));
}

//-------------------------
void ParticleFactory::createChannel(DecayingParticle* parent, Particle* daughterA, Particle* daughterB, unsigned L)
{
    yap::SpinAmplitude* amplitude = new yap::SpinAmplitude(parent->quantumNumbers(), daughterA->quantumNumbers(), daughterB->quantumNumbers());
    parent->addChannel(new yap::DecayChannel(daughterA, daughterB, L, *amplitude));
}

//-------------------------
QuantumNumbers ParticleFactory::createQuantumNumbers(int PDG)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);

    unsigned char twoJ = std::round(p->Spin() * 2.);
    signed char P(p->Parity());

    // C-parity
    signed char C(0);
    if (p->Charge() == 0) { // only defined for neutral states
        if (p == p->AntiParticle())
            C = 1;
        else
            C = -1;
    }

    signed char I = p->Isospin();

    // TODO: G-Parity
    signed char G = 0;

    //std::cout << "p->Spin() " << p->Spin() << " twoJ " << (int)twoJ << " C " << (int)C << "\n";
    //std::cout << QuantumNumbers(twoJ, P, C, I, G);

    return QuantumNumbers(twoJ, P, C, I, G);
}

}
