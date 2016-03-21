#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <MassAxes.h>
#include <Model.h>
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
 * Downstream theta angles must be the same regardless of rotations/boosts of the final state particles
 */

TEST_CASE( "HelicityAngles_boostRotate" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
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
    rho->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

    // create DataSet
    auto data = M.dataSet();

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    std::vector<double> masses = { piPlus->mass()->value(), piMinus->mass()->value(), piPlus->mass()->value() };


    for (unsigned int iEvt = 0; iEvt < 1000; ++iEvt) {
        std::map<std::shared_ptr<yap::ParticleCombination>, std::vector<double>> resultingThetas;

        TGenPhaseSpace event;
        event.SetDecay(P, masses.size(), &masses[0]);
        event.Generate();

        for (unsigned iTrans=0; iTrans<7; ++iTrans) {

            std::vector<yap::FourVector<double> > momenta;
            std::vector<TLorentzVector> root_momenta;

            double angle = gRandom->Uniform(0, yap::pi<double>());
            double boost = gRandom->Uniform(-0.99, 0.99);

            for (unsigned i = 0; i < masses.size(); ++i) {
                TLorentzVector p = *event.GetDecay(i);

                // testing. Theta of downstream helicity angles must stay the same
                // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes

                switch (iTrans) {
                case 0:
                default:
                    break;
                case 1:
                    //std::cout << "rotate around X by " << angle << "\n";
                    p.RotateX(angle);
                    break;
                case 2:
                    //std::cout << "rotate around Y by " << angle << "\n";
                    p.RotateY(angle);
                    break;
                case 3:
                    //std::cout << "rotate around Z by " << angle << "\n";
                    p.RotateZ(angle);
                    break;
                case 4:
                    //std::cout << "boost along X by " << boost << "\n";
                    p.Boost(boost, 0., 0.);
                    break;
                case 5:
                    //std::cout << "boost along Y by " << boost << "\n";
                    p.Boost(0., boost, 0.);
                    break;
                case 6:
                    //std::cout << "boost along Z by " << boost << "\n";
                    p.Boost(0., 0., boost);
                    break;
                }

                momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));
                root_momenta.push_back(p);
            }

            data.add(momenta);
            const auto dp = data.points().back();

            // check 4 momenta in their rest frames
            /*for (auto& pc : M.helicityAngles()->particleCombinations()) {
                auto root_momenta_copy(root_momenta);
                TLorentzVector b;
                for (unsigned i : pc->indices())
                    b += root_momenta_copy[i];

                std::cout << "4momenta in " << to_string(*pc) << " rest frame:\n";
                for (auto& v : root_momenta_copy) {
                    v.Boost(-b.BoostVector());
                    v.Print();
                }

            }*/

            // compare results
            //std::cout << "transformation " << iTrans << "\n";
            for (auto& pc : M.helicityAngles()->particleCombinations()) {
                //std::cout << std::left << std::setw(50) << yap::to_string(*pc);
                //std::cout << "  helicityAngles():  (" <<  M.helicityAngles()->phi(dp, pc) << ", " << M.helicityAngles()->theta(dp, pc) << ")\n";

                if (pc->indices().size() < M.finalStateParticles().size())
                    resultingThetas[pc].push_back(M.helicityAngles()->theta(dp, pc));
            }

            //std::cout << "\n";
        }


        // check if thetas are equal
        // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
        for (auto& kv : resultingThetas) {
            double theta_0 = kv.second[0];
            for (auto theta : kv.second) {
                //std::cout << to_string(*kv.first) << ": " << theta << "\n";
                REQUIRE(theta == Approx(theta_0));
            }
        }

        //std::cout << "-------------------------------------------------------\n";
    }
}
