#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <Attributes.h>
#include <logging.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PHSP.h>
#include <Rotation.h>

#include <assert.h>
#include <cmath>
#include <random>

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

    yap::ParticleTable T = yap::read_pdl_file(find_pdl_file());

    double D_mass = T["D+"].mass();

    auto M = d3pi<yap::HelicityFormalism>();

    // choose default Dalitz coordinates
    auto A = M->massAxes();
    // get mass^2 ranges
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M->finalStateParticles()));

    // create DataSet
    auto data = M->createDataSet();

    // create random number engine for generation of points
    std::mt19937 g(0);

    // create random number generators
    std::uniform_real_distribution<double> uniform_0_pi(0, yap::pi());
    std::uniform_real_distribution<double> uniform_m99_p99(-0.99, 0.99);

    for (unsigned int iEvt = 0; iEvt < 1000; ++iEvt) {
        yap::ParticleCombinationMap<std::vector<double> > resultingThetas;

        // generate random phase space point (with 100 attempts before failing)
        auto momenta = yap::phsp(*M, D_mass, A, m2r, g, 100);
        if (momenta.empty())
            continue;

        for (unsigned iTrans = 0; iTrans < 7; ++iTrans) {

            double angle = uniform_0_pi(g);
            double boost = uniform_m99_p99(g);

            for (auto& p : momenta) {

                // testing. Theta of downstream helicity angles must stay the same
                switch (iTrans) {
                    // rotate around axis: case 1, x; case 2, y; case 3, z
                    case 1:
                        p = yap::rotation(yap::ThreeVector<double>({1., 0., 0.}), angle) * p;
                        break;
                    case 2:
                        p = yap::rotation(yap::ThreeVector<double>({0., 1., 0.}), angle) * p;
                        break;
                    case 3:
                        p = yap::rotation(yap::ThreeVector<double>({0., 0., 1.}), angle) * p;
                        break;
                    // boost in direction of axis: case 4, x; case 5, y; case 6, z
                    case 4:
                        p = yap::lorentzTransformation(yap::ThreeVector<double>({boost, 0, 0})) * p;
                        break;
                    case 5:
                        p = yap::lorentzTransformation(yap::ThreeVector<double>({0, boost, 0})) * p;
                        break;
                    case 6:
                        p = yap::lorentzTransformation(yap::ThreeVector<double>({0, 0, boost})) * p;
                        break;
                    default:
                        break;
                }
            }

            data.push_back(momenta);
            const auto dp = data.back();

            auto rho = particle(*M, yap::is_named("rho0"));
            
            // compare results
            for (auto& pc_rho : rho->particleCombinations())
                resultingThetas[pc_rho].push_back(M->helicityAngles()(dp, data, pc_rho).theta);
        }

        // check if thetas are equal
        // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
        for (auto& kv : resultingThetas) {
            double theta_0 = kv.second[0];
            for (auto theta : kv.second)
                REQUIRE(theta == Approx(theta_0).epsilon(0.001));
        }
    }
}
