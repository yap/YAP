#include <catch.hpp>
#include <catch_capprox.hpp>

#include <BreitWigner.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <ZemachSpinAmplitude.h>

#include <cmath>

TEST_CASE( "swapDalitzAxes" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
    //yap::plainLogs(el::Level::Global);


    std::vector<double> m2({1.4, 1.5});


    std::vector<std::vector<double> > parameters{m2, {m2[1], m2[0]}};

    std::vector<double> resultingAmplitudes({0,0});

    for (unsigned i=0; i<2; ++i) {

        auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

        // final state particles
        auto kPlus  = F.fsp(+321);
        auto kMinus = F.fsp(-321);
        auto piPlus = F.fsp(+211);

        auto M = std::make_unique<yap::Model>(std::make_unique<yap::ZemachFormalism>());

        //M->setFinalState({kPlus, kMinus, piPlus});
        M->setFinalState({piPlus, kMinus, kPlus});

        // use common radial size for all resonances
        double radialSize = 3.; // [GeV^-1]

        // initial state particle
        auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

        auto piK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.75, "piK0", radialSize, std::make_shared<yap::BreitWigner>(0.75, 0.025));
        piK0->addChannel({piPlus, kMinus});
        D->addChannel({piK0, kPlus})->freeAmplitudes()[0]->setValue(0.5 * yap::Complex_1);

        auto piK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.00, "piK1", radialSize, std::make_shared<yap::BreitWigner>(1.00, 0.025));
        piK1->addChannel({piPlus, kMinus});
        D->addChannel({piK1, kPlus})->freeAmplitudes()[0]->setValue(1. * yap::Complex_1);

        auto piK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.25, "piK2", radialSize, std::make_shared<yap::BreitWigner>(1.25, 0.025));
        piK2->addChannel({piPlus, kMinus});
        D->addChannel({piK2, kPlus})->freeAmplitudes()[0]->setValue(30. * yap::Complex_1);

        M->initializeForMonteCarloGeneration(1);

        // choose Dalitz coordinates m^2_12 and m^2_23
        yap::MassAxes massAxes;
        if (i==0)
            massAxes = M->getMassAxes({{0, 1}, {1, 2}});
        else
            massAxes = M->getMassAxes({{2, 1}, {1, 0}});


        // calculate four-momenta
        auto P = M->calculateFourMomenta(massAxes, parameters[i]);

        // if failed, outside phase space
        if (P.empty()) {
            WARN( "outside of phasespace" );
            break;
        }

        M->setFinalStateMomenta(M->dataSet()[0], P, 0);

        resultingAmplitudes[i] = M->logOfSquaredAmplitude(M->dataSet()[0], 0);

    }

    REQUIRE( resultingAmplitudes[0] == resultingAmplitudes[1]);




}
