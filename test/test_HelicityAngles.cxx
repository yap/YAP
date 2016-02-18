#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicitySpinAmplitude.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <MathUtilities.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <WignerD.h>
#include <ZemachSpinAmplitude.h>

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

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // set final state
    M.setFinalState({piPlus, piMinus, piPlus});

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.getMassAxes({{0, 1}, {1, 2}});

    REQUIRE( M.consistent() );

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

        M.addDataPoint(momenta);
        auto dp = M.dataSet().back();


        // \todo we are now assuming that the D's coordinate system is (1,0,0; 0,1,0; 0,0,1)
        // this is ok since it is the default and we didn't change it
        // actually we would have to rotate

        for (auto pc : D->particleCombinations()) {
            REQUIRE( M.fourMomenta().m(dp, pc) == Approx(D->mass()->value()) );

            TLorentzVector daughter;
            for (unsigned i : pc->daughters()[0]->indices()) {
                daughter += *event.GetDecay(i);
            }

            DEBUG( to_string(*pc) << "; phi = " << daughter.Phi() << "; theta = " << daughter.Theta());

            REQUIRE( M.helicityAngles().phi(dp, pc) == Approx(daughter.Phi()) );
            REQUIRE( M.helicityAngles().theta(dp, pc) == Approx(daughter.Theta()) );

            // \todo recurse down
        }

    }
}
