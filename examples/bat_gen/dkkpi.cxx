// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "dkkpi.h"

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include "ParticleCombination.h"
#include "QuantumNumbers.h"
#include "Resonance.h"

#include <complex>

std::unique_ptr<yap::Model> dkkpi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{

    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // final state particles
    auto kPlus  = F.fsp(+321);
    auto kMinus = F.fsp(-321);
    auto piPlus = F.fsp(+211);

    auto M = std::make_unique<yap::Model>(std::move(SAC));

    M->setFinalState({kPlus, kMinus, piPlus});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    auto KK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 1.1, "KK0", radialSize, std::make_shared<yap::BreitWigner>(1.1, 0.05));
    KK0->addChannel({kPlus, kMinus});
    D->addChannel({KK0, piPlus})->freeAmplitudes()[0]->setValue(yap::Complex_1);

    auto KK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.35, "KK1", radialSize, std::make_shared<yap::BreitWigner>(1.35, 0.05));
    KK1->addChannel({kPlus, kMinus});
    D->addChannel({KK1, piPlus})->freeAmplitudes()[0]->setValue(2. * yap::Complex_1);

    auto KK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.6, "KK2", radialSize, std::make_shared<yap::BreitWigner>(1.6, 0.05));
    KK2->addChannel({kPlus, kMinus});
    D->addChannel({KK2, piPlus})->freeAmplitudes()[0]->setValue(30. * yap::Complex_1);

    /*
    // X_2
    auto X_2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.2, "X_2", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(X_2->massShape())->width()->setValue(80e-3);
    X_2->addChannel({piPlus, kMinus});
    D->addChannel({X_2, kPlus});
    */

    return M;
}
