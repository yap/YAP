#include <catch.hpp>
#include <catch_capprox.hpp>

#include "BreitWigner.h"
#include "Constants.h"
#include "DataPoint.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "HelicityFormalism.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "Parameter.h"
#include "Particle.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PHSP.h"
#include "Resonance.h"
#include "Spin.h"
#include "UnitSpinAmplitude.h"

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

    /// Constructor
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    testHelicitySpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
        HelicitySpinAmplitude(two_J, two_j, l, two_s)
    { }
};

}

std::shared_ptr<yap::Model> create_model()
{
    auto M = std::make_shared<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    double radialSize = 1.;

    // initial state particle
    auto D = factory.decayingParticle(421, radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // Set final-state particles
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = factory.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = factory.resonance(20213, radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addChannel(sigma, piPlus);
    a_1->addChannel(rho,   piPlus);

    // D's channels
    D->addChannel(rho, rho);
    D->addChannel(a_1, piMinus);

    M->addInitialStateParticle(D);

    return M;
}

yap::DataSet generate_data(std::shared_ptr<yap::Model> M)
{
    auto isp = M->initialStateParticles().begin()->first;
    auto A = M->massAxes();
    auto m2r = yap::squared(yap::mass_range(A, isp, M->finalStateParticles()));

    yap::DataSet data(M->createDataSet());

    std::mt19937 g(0);
    // fill data set with 1 point
    std::generate_n(std::back_inserter(data), 1,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*M), isp->mass()->value(), A, m2r, g, std::numeric_limits<unsigned>::max()));

    return data;
}

TEST_CASE( "UnitSpinAmplitude" )
{
    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto M = create_model();
    auto data = generate_data(M);
    auto isp = M->initialStateParticles().begin()->first;

    yap::HelicityFormalism formalism;

    static const unsigned J_max = 4;
    unsigned false_unit(0);

    for (unsigned two_J = 0; two_J < 2*J_max; two_J += 2)
        for (unsigned two_j2 = 0; two_j2 < 2*J_max; two_j2 += 2)
            for (unsigned two_j1 = 0; two_j1 < 2*J_max; two_j1 += 2)
                for (unsigned two_s = 0; two_s < 2*J_max; two_s += 2)
                    for (unsigned l = 0; l < J_max; ++l) {

                        const yap::SpinVector two_j({two_j1, two_j2});

                        try {
                            auto sa = formalism.spinAmplitude(two_J, two_j, l, two_s);
                            auto sa_test = yap::testHelicitySpinAmplitude(two_J, two_j, l, two_s);
                            sa_test.setModel(*M);

                            DEBUG("two_J = " << two_J << "; two_j1 = " << two_j1 << "; two_j2 = " << two_j2
                                    << "; two_s = " << two_s << "; l = " << l);

                            for (auto pc_cache : M->particleCombinationCache()) {
                                auto pc = pc_cache.lock();
                                if (not pc
                                        or pc->isFinalStateParticle()
                                        or not is_initial_state_particle_combination(*pc->origin(), M.get()))
                                    continue;

                                //DEBUG("pc " << to_string_with_parent(*pc));

                                REQUIRE( sa->twoM() == sa_test.twoM() );

                                // loop over spin projections
                                for (int two_M : sa_test.twoM()) {

                                    REQUIRE( sa->twoM(two_M) == sa_test.twoM(two_M) );

                                    for (auto two_m : sa_test.twoM(two_M)) {
                                        auto amp = sa_test.calc(two_M, two_m, data[0], pc);
                                        if (std::dynamic_pointer_cast<yap::UnitSpinAmplitude>(sa))
                                            // if it is a UnitSpinAmplitude, check if the corresponding HelicitySpinAmplitude is 1
                                            REQUIRE(amp == Catch::Detail::CApprox(yap::Complex_1));
                                        else {
                                            // if it is NOT a UnitSpinAmplitude, check if the corresponding HelicitySpinAmplitude is != 1
                                            //REQUIRE(amp != yap::Complex_1);
                                            if (amp == yap::Complex_1 and sa->size() > 0) {
                                                DEBUG("is 1 but not Unit and has size > 0");
                                                ++false_unit;
                                            }
                                        }

                                    }
                                }
                            }
                        }
                        catch (const yap::exceptions::AngularMomentumNotConserved&) { /* ignore */ }
                    }

    // if false_unit is a big number, one can probably make UnitSpinAmplitudes in more cases
    DEBUG("In " << false_unit << " cases the amplitude was 1, but it was not a UnitSpinAmplitude.");
    REQUIRE(false_unit < 3);

}

