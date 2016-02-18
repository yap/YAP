#include <catch.hpp>
#include <catch_capprox.hpp>

#include <BreitWigner.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <ZemachFormalism.h>

#include <cmath>

/**
 * Test that the amplitude remains the same after swapping the order of the final state particles
 */

TEST_CASE( "swapFinalStates" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
    //yap::plainLogs(el::Level::Global);

    std::array<double, 4> PDGs = {411, 321, -321, 211}; // D+ -> K+ K- pi+
    std::array<double, 2> m2_ab_range = {0.4, 1.9};
    std::array<double, 2> m2_bc_range = {0.9, 3.1};

    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");


    const unsigned N = 20;
    for (double m2_ab = m2_ab_range[0]; m2_ab <= m2_ab_range[1]; m2_ab += (m2_ab_range[1] - m2_ab_range[0]) / N) {
        for (double m2_bc = m2_bc_range[0]; m2_bc <= m2_bc_range[1]; m2_bc += (m2_bc_range[1] - m2_bc_range[0]) / N) {

            // calc 3rd inv mass square
            double m2_ac = pow(F.decayingParticle(PDGs[0], 3.)->mass()->value(), 2)
                           + pow(F.fsp(PDGs[1])->mass()->value(), 2)
                           + pow(F.fsp(PDGs[2])->mass()->value(), 2)
                           + pow(F.fsp(PDGs[3])->mass()->value(), 2)
                           - m2_ab - m2_bc;

            if (m2_ac < 0.) {
                //std::cout << "m2_ac < 0.\n";
                continue;
            }


            // loop over SpinFormalisms
            for (unsigned i_formalism = 0; i_formalism < 2; ++i_formalism) {

                std::vector<double> resultingAmplitudes(6, 0.);

                std::cout << (i_formalism == 0 ? "1. Zemach  \t" : "2. Helicity\t") << "    ";

                // loop over fsp swaps
                for (unsigned i = 0; i < 6; ++i) {

                    // final state particles
                    auto kPlus  = F.fsp(PDGs[1]);
                    auto kMinus = F.fsp(PDGs[2]);
                    auto piPlus = F.fsp(PDGs[3]);

                    auto M = i_formalism == 0 ? std::make_unique<yap::Model>(std::make_unique<yap::ZemachFormalism>()) :
                             std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());


                    switch (i) {
                        case 0:
                            // original
                            M->setFinalState({piPlus, kMinus, kPlus});
                            break;
                        case 1:
                            // 0 <-> 1
                            M->setFinalState({kMinus, piPlus, kPlus});
                            break;
                        case 2:
                            // 0 <-> 2
                            M->setFinalState({kPlus, kMinus, piPlus});
                            break;
                        case 3:
                            // 1 <-> 2
                            M->setFinalState({piPlus, kPlus, kMinus});
                            break;
                        case 4:
                            // <<
                            M->setFinalState({kMinus, kPlus, piPlus});
                            break;
                        case 5:
                            // >>
                            M->setFinalState({kPlus, piPlus, kMinus});
                            break;
                    }


                    // initial state particle
                    auto D = F.decayingParticle(PDGs[0], 3.);

                    /*auto piK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.75, "piK0", 3., std::make_shared<yap::BreitWigner>(0.75, 0.025));
                    piK0->addChannel({piPlus, kMinus});
                    D->addChannel({piK0, kPlus})->freeAmplitudes()[0]->setValue(0.5 * yap::Complex_1);
                    */
                    auto piK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.00, "piK1", 3., std::make_shared<yap::BreitWigner>(1.00, 0.025));
                    piK1->addChannel({piPlus, kMinus});
                    D->addChannel({piK1, kPlus})->freeAmplitudes()[0]->setValue(1. * yap::Complex_1);
                    /*
                                        auto piK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.25, "piK2", 3., std::make_shared<yap::BreitWigner>(1.25, 0.025));
                                        piK2->addChannel({piPlus, kMinus});
                                        D->addChannel({piK2, kPlus})->freeAmplitudes()[0]->setValue(30. * yap::Complex_1);
                    */
                    M->initializeForMonteCarloGeneration(1);


                    // Dalitz coordinates
                    yap::MassAxes massAxes = M->getMassAxes({{0, 1}, {1, 2}});;
                    std::vector<double> squared_masses;
                    switch (i) {
                        case 0:
                            // original
                            squared_masses = {m2_ab, m2_bc};
                            break;
                        case 1:
                            // 0 <-> 1
                            squared_masses = {m2_ab, m2_ac};
                            break;
                        case 2:
                            // 0 <-> 2
                            squared_masses = {m2_bc, m2_ab};
                            break;
                        case 3:
                            // 1 <-> 2
                            squared_masses = {m2_ac, m2_bc};
                            break;
                        case 4:
                            // <<
                            squared_masses = {m2_bc, m2_ac};
                            break;
                        case 5:
                            // >>
                            squared_masses = {m2_ac, m2_ab};
                            break;
                    }

                    // calculate four-momenta
                    auto P = M->calculateFourMomenta(massAxes, squared_masses);

                    // if failed, outside phase space
                    if (P.empty()) {
                        std::cout << "PhSp    ";
                        continue;
                    }

                    /// \todo enabling this makes a difference in the amplitude, but it should not
                    //P = yap::operator-(lorentzTransformation(P)) * P;
                    M->setFinalStateMomenta(M->dataSet()[0], P, 0);

                    resultingAmplitudes[i] = M->logOfSquaredAmplitude(M->dataSet()[0], 0);

                    std::cout << resultingAmplitudes[i] << "   ";
                }

                /*REQUIRE( std::isfinite(resultingAmplitudes[0]) );
                for (unsigned i=1; i<resultingAmplitudes.size(); ++i)
                    REQUIRE( resultingAmplitudes[0] == Approx(resultingAmplitudes[i]) );*/

                std::cout << "\n";

            }// end loop over SpinFormalisms

        }
    }




}