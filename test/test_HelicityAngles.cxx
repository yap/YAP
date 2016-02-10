#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <HelicitySpinAmplitude.h>
#include <InitialStateParticle.h>
#include <logging.h>
#include <MassAxes.h>
#include <MathUtilities.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <WignerD.h>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <cmath>

TEST_CASE( "HelicityAngles" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
    //yap::plainLogs(el::Level::Debug);

    // init random generator
    gRandom->SetSeed(1234);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

    // initial state particle
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize,
            //std::make_unique<yap::ZemachSpinAmplitudeCache>());
            std::make_unique<yap::HelicitySpinAmplitudeCache>());

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(211);
    std::shared_ptr<yap::FinalStateParticle> piMinus = factory.createFinalStateParticle(-211);

    // set final state
    D->setFinalStateParticles({piPlus, piMinus, piPlus});

    // rho
    std::shared_ptr<yap::Resonance> rho = std::make_shared<yap::Resonance>(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(rho->massShape()).width()->setValue(0.149);
    rho->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = D->getMassAxes({{0, 1}, {1, 2}});

    D->prepare();

    const int nParticles(3);

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[nParticles] = { piPlus->mass()->value(), piMinus->mass()->value(),
                                    piPlus->mass()->value()
                                  };

    for (unsigned int iEvt = 0; iEvt < 10; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, nParticles, masses);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        for (unsigned i = 0; i < nParticles; ++i) {
            TLorentzVector p = *event.GetDecay(i);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));
        }

        D->addDataPoint(momenta);
        auto dp = D->dataSet().back();


        // \todo we are now assuming that the D's coordinate system is (1,0,0; 0,1,0; 0,0,1)
        // actually we would have to rotate

        for (auto pc : D->particleCombinations()) {
            TLorentzVector daughter;
            for (unsigned i : pc->daughters()[0]->indices())
                daughter += *event.GetDecay(i);

            DEBUG( to_string(*pc) << "; phi = " << daughter.Phi() << "; theta = " << daughter.Theta());

            REQUIRE( D->helicityAngles().phi(dp, pc) == Approx(daughter.Phi()) );
            REQUIRE( D->helicityAngles().theta(dp, pc) == Approx(daughter.Theta()) );

            // \todo recurse down
        }

    }
}
