#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>

#include <TGenPhaseSpace.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <assert.h>
#include <cmath>

/*
 * Test the calculation of helicity angles
 * Downstream theta angles must be the same regardless of rotations/boosts of the final state particles.
 * Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
 */

TEST_CASE( "HelicityAngles_boostRotate" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    // init random generator
    gRandom->SetSeed(1234);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

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
    rho->addChannelSA({piPlus, piMinus});

    // Add channels to D
    D->addChannelSA({rho,      piPlus});

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

    // create DataSet
    auto data = M.createDataSet();

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    std::vector<double> masses = { piPlus->mass()->value(), piMinus->mass()->value(), piPlus->mass()->value() };

    for (unsigned int iEvt = 0; iEvt < 1000; ++iEvt) {
        std::map<std::shared_ptr<yap::ParticleCombination>, std::vector<double>> resultingThetas;

        TGenPhaseSpace event;
        event.SetDecay(P, masses.size(), &masses[0]);
        event.Generate();

        for (unsigned iTrans = 0; iTrans < 7; ++iTrans) {

            std::vector<yap::FourVector<double> > momenta;
            std::vector<TLorentzVector> root_momenta;

            double angle = gRandom->Uniform(0, yap::pi<double>());
            double boost = gRandom->Uniform(-0.99, 0.99);

            for (unsigned i = 0; i < masses.size(); ++i) {
                TLorentzVector p = *event.GetDecay(i);

                // testing. Theta of downstream helicity angles must stay the same
                switch (iTrans) {
                    case 0:
                    default:
                        break;
                    case 1:
                        p.RotateX(angle);
                        break;
                    case 2:
                        p.RotateY(angle);
                        break;
                    case 3:
                        p.RotateZ(angle);
                        break;
                    case 4:
                        p.Boost(boost, 0., 0.);
                        break;
                    case 5:
                        p.Boost(0., boost, 0.);
                        break;
                    case 6:
                        p.Boost(0., 0., boost);
                        break;
                }

                momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));
                root_momenta.push_back(p);
            }

            data.add(momenta);
            const auto dp = data.points().back();

            // compare results
            for (auto& pc : M.helicityAngles()->particleCombinations())
                if (pc->indices().size() < M.finalStateParticles().size())
                    resultingThetas[pc].push_back(M.helicityAngles()->theta(dp, pc));
        }

        // check if thetas are equal
        // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
        for (auto& kv : resultingThetas) {
            double theta_0 = kv.second[0];
            for (auto theta : kv.second)
                REQUIRE(theta == Approx(theta_0));
        }
    }
}
