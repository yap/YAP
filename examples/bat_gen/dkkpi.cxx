// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "dkkpi.h"

#include <BAT/BCMath.h>

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "HelicitySpinAmplitude.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "QuantumNumbers.h"
#include "Resonance.h"
#include "ZemachSpinAmplitude.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <complex>
#include <iostream>

// ---------------------------------------------------------
dkkpi::dkkpi(std::string name)
    : d2X(name)
{
    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> kPlus  = Factory_.createFinalStateParticle(+321);
    std::shared_ptr<yap::FinalStateParticle> kMinus = Factory_.createFinalStateParticle(-321);
    std::shared_ptr<yap::FinalStateParticle> piPlus = Factory_.createFinalStateParticle(+211);

    // initial state particle
    D_ = Factory_.createInitialStateParticle(Factory_.pdgCode("D+"), radialSize, std::make_unique<yap::ZemachSpinAmplitudeCache>());
    // D_ = Factory_.createInitialStateParticle(Factory_.pdgCode("D+"), radialSize, std::make_unique<yap::HelicitySpinAmplitudeCache>());
    D_->setFinalStateParticles({kPlus, kMinus, piPlus});

    // phi
    // std::shared_ptr<yap::Resonance> phi = std::make_shared<yap::Resonance>(Factory_.quantumNumbers("phi"), 1010.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    std::shared_ptr<yap::Resonance> phi = std::make_shared<yap::Resonance>(yap::QuantumNumbers(2, 0), 1310.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(phi->massShape()).width()->setValue(20e-3);
    phi->addChannel({kPlus, kMinus});
    D_->addChannel({phi, piPlus});

    std::cout << *phi << std::endl;

    /*
    // X_2
    std::shared_ptr<yap::Resonance> X_2 = std::make_shared<yap::Resonance>(yap::QuantumNumbers(4, 0), 1.2, "X_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(X_2->massShape()).width()->setValue(80e-3);
    X_2->addChannel({piPlus, kMinus});
    D_->addChannel({X_2, kPlus});
    */

    D_->prepare();

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

}
