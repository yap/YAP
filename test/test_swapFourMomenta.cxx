#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
#include <logging.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <PDL.h>

#include <cmath>
#include <assert.h>

/**
 * Test that the amplitude remains the same after swapping the four momenta of the final state particles
 */

TEST_CASE( "swapFourMomenta" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto m = d4pi();
    const double D0_mass = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl")["D0"].mass();

    const auto massAxes = m->massAxes();
    auto m2r = yap::squared(mass_range(D0_mass, m->massAxes(), m->finalStateParticles()));
    const auto FSPs = m->finalStateParticles();
    const auto components = m->components();
    LOG(INFO) << components.size();
    assert(components.size() == 1);

    std::mt19937 g(164852419);

    for (unsigned i = 0; i < 100; ++i) {
        std::vector<yap::FourVector<double> > P = yap::phsp(*m, D0_mass, massAxes, m2r, g, std::numeric_limits<unsigned>::max());

        yap::DataSet dataSet(*m);

        dataSet.push_back({P[0], P[1], P[2], P[3]});
        dataSet.push_back({P[2], P[1], P[0], P[3]}); // 0 <-> 2
        dataSet.push_back({P[0], P[3], P[2], P[1]}); // 1 <-> 2
        dataSet.push_back({P[2], P[3], P[0], P[1]}); // both

        m->calculate(dataSet);

        std::vector<std::complex<double> > amplitudes;

        for (auto& dp : dataSet)
            amplitudes.push_back(yap::amplitude(components[0].decayTrees(), dp));

        for (unsigned i = 1; i < amplitudes.size(); ++i) 
            REQUIRE(amplitudes[i] == Catch::Detail::CApprox(amplitudes[0]));
    }

}
