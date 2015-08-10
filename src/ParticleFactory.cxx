#include "ParticleFactory.h"
#include "BreitWigner.h"
#include "FinalStateParticle.h"
#include "Resonance.h"

#include <TDatabasePDG.h>
#include "logging.h"

namespace yap {

//-------------------------
Particle* ParticleFactory::createFinalStateParticle(int PDG)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new FinalStateParticle(createQuantumNumbers(PDG), p->Mass(), std::string(p->GetName()), PDG);
}

//-------------------------
Particle* ParticleFactory::createResonance(int PDG, double radialSize, const MassShape& massShape)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new Resonance(createQuantumNumbers(PDG), p->Mass(), p->GetName(), radialSize, massShape);
}

//-------------------------
Particle* ParticleFactory::createResonanceBreitWigner(int PDG, double radialSize)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new Resonance(createQuantumNumbers(PDG), p->Mass(), std::string(p->GetName()),
        radialSize, BreitWigner(p->Mass(), p->Width()));
}

//-------------------------
QuantumNumbers ParticleFactory::createQuantumNumbers(int PDG)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);

    unsigned char twoJ = std::round(p->Spin() * 2.);

    // C-parity
    char C = 0;
    if (p->Charge() != 0) { // only defined for neutral states
        if (p == p->AntiParticle())
            C = 1;
        else
            C = -1;
    }

    // TODO: G-Parity
    char G = 0;

    return QuantumNumbers(twoJ, p->Parity(), C, p->Isospin(), G);
}

}
