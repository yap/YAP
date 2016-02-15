#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
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

#include <cmath>

TEST_CASE( "FourMomentaCalculation" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
    //yap::plainLogs(el::Level::Global);


    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

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

    M.initializeForMonteCarloGeneration(1);
    REQUIRE( M.consistent() );

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.getMassAxes({{0, 1}, {1, 2}});

    std::vector<double> m2 = {0.9, 1.1};

    // calculate final state momenta
    auto P = M.calculateFourMomenta(massAxes, m2);


    //
    // require invariant masses to be correct
    //
    // fsp
    for (auto& v : P)
        REQUIRE( abs(v) == Approx(piPlus->mass()->value()) );
    // isp
    yap::FourVector<double> isp = std::accumulate(P.begin(), P.end(), yap::FourVector_0);
    REQUIRE( abs(isp) == Approx(D->mass()->value()) );
    // 2-p
    unsigned i_a(0);
    for (auto& axis : massAxes) {
        auto sum = yap::FourVector_0;
        for (unsigned i : axis->indices())
            sum +=  P[i];

        REQUIRE( abs(sum) == Approx(m2[i_a]) );
        ++i_a;
    }


    // Require ISP to be in its rest frame
    REQUIRE( abs(vect(isp)) == Approx(0) );
}
