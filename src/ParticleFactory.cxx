#include "ParticleFactory.h"
#include "FinalStateParticle.h"

#include <TDatabasePDG.h>
#include "logging.h"

namespace yap {

//-------------------------
Particle* ParticleFactory::createFinalStateParticle(int PDG)
{
    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(PDG);
    return new FinalStateParticle(createQuantumNumbers(PDG), p->Mass(), p->GetName(), PDG);
}

//-------------------------
Particle* ParticleFactory::createResonance(int PDG)
{

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
