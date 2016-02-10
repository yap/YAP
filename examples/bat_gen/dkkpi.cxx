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
    auto kPlus  = Factory_.fsp(+321);
    auto kMinus = Factory_.fsp(-321);
    auto piPlus = Factory_.fsp(+211);

    // initial state particle
    D_ = Factory_.isp(Factory_.pdgCode("D+"), radialSize, std::make_unique<yap::ZemachSpinAmplitudeCache>());
    // D_ = Factory_.isp(Factory_.pdgCode("D+"), radialSize, std::make_unique<yap::HelicitySpinAmplitudeCache>());
    D_->setFinalState({kPlus, kMinus, piPlus});

    // phi
    // auto phi = std::make_shared<yap::Resonance>(Factory_.quantumNumbers("phi"), 1010.e-3, "phi", radialSize, std::make_shared<yap::BreitWigner>());
    auto phi = std::make_shared<yap::Resonance>(yap::QuantumNumbers(2, 0), 1310.e-3, "phi", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(phi->massShape())->width()->setValue(20e-3);
    phi->addChannel({kPlus, kMinus});
    D_->addChannel({phi, piPlus});

    std::cout << *phi << std::endl;

    /*
    // X_2
    auto X_2 = std::make_shared<yap::Resonance>(yap::QuantumNumbers(4, 0), 1.2, "X_2", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(X_2->massShape())->width()->setValue(80e-3);
    X_2->addChannel({piPlus, kMinus});
    D_->addChannel({X_2, kPlus});
    */

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

}
