#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <DataPoint.h>
#include <logging.h>
#include <MathUtilities.h>
#include <Particle.h>
#include <Spin.h>
#include <UnitSpinAmplitude.h>

#include <assert.h>
#include <memory>
#include <random>

/**
 *  Test that when HelicityFormalism decides to use a UnitSpinAmplitude, that
 *  the HelicitySpinAmplitude with the same quantum numbers also calculates 1
 */


namespace yap {

// to access the protected constructor publicly
class testHelicitySpinAmplitude : public HelicitySpinAmplitude
{
public:

    testHelicitySpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
        HelicitySpinAmplitude(m, two_J, two_j, l, two_s)
    { }
};

}

//-------------------------
TEST_CASE( "UnitSpinAmplitude" )
{
    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto M = create_model();
    auto data = generate_data(M, 1);
    auto isp = M->initialStateParticles().begin()->first;

    yap::HelicityFormalism formalism;

    static const unsigned J_max = 4;
    unsigned false_unit(0);

    for (unsigned two_J = 0; two_J < 2 * J_max; two_J += 2)
        for (unsigned two_j2 = 0; two_j2 < 2 * J_max; two_j2 += 2)
            for (unsigned two_j1 = 0; two_j1 < 2 * J_max; two_j1 += 2)
                for (unsigned two_s = 0; two_s < 2 * J_max; two_s += 2)
                    for (unsigned l = 0; l < J_max; ++l) {

                        const yap::SpinVector two_j({two_j1, two_j2});

                        try {
                            auto sa = formalism.spinAmplitude(two_J, two_j, l, two_s);
                            auto sa_test = yap::testHelicitySpinAmplitude(*M, two_J, two_j, l, two_s);

                            DEBUG("two_J = " << two_J << "; two_j1 = " << two_j1 << "; two_j2 = " << two_j2
                                  << "; two_s = " << two_s << "; l = " << l);

                            REQUIRE( sa->twoM() == sa_test.twoM() );

                            for (auto two_M : sa_test.twoM())
                                REQUIRE( sa->twoM(two_M) == sa_test.twoM(two_M) );

                            for (auto pc_cache : M->particleCombinationCache()) {
                                auto pc = pc_cache.lock();
                                if (not pc
                                        or pc->isFinalStateParticle()
                                        or not is_initial_state_particle_combination(*pc->origin(), M.get()))
                                    continue;

                                // loop over spin projections
                                for (int two_M : sa_test.twoM()) {

                                    for (auto two_m : sa_test.twoM(two_M)) {

                                        auto amp = sa_test.calc(two_M, two_m, data[0], pc);

                                        if (std::dynamic_pointer_cast<yap::UnitSpinAmplitude>(sa))
                                            // if it is a UnitSpinAmplitude, check if the corresponding HelicitySpinAmplitude is 1
                                            REQUIRE(amp == Catch::Detail::CApprox(1.));

                                        else {
                                            // if it is NOT a UnitSpinAmplitude, check if the corresponding HelicitySpinAmplitude is != 1
                                            if (amp == 1 and sa->size() > 0) {
                                                DEBUG("is 1 but not Unit and has size > 0");
                                                ++false_unit;
                                            }
                                        }

                                    }
                                }
                            }
                        } catch (const yap::exceptions::AngularMomentumNotConserved&) { /* ignore */ }
                    }

    // if false_unit is a big number, one can probably make UnitSpinAmplitudes in more cases
    DEBUG("In " << false_unit << " cases the amplitude was 1, but it was not a UnitSpinAmplitude.");
    REQUIRE(false_unit < 3);

}

