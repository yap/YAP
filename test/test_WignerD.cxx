#include <catch.hpp>
#include <catch_capprox.hpp>

#include <Exceptions.h>
#include <logging.h>
#include <MathUtilities.h>
#include <WignerD.h>

#include <cmath>

void checkDSymmetries(unsigned twoJ, int twoM, int twoN, double alpha, double beta, double gamma)
{
    REQUIRE( yap::DFunction(twoJ, twoM, twoN, alpha, beta, gamma) == std::conj(yap::DFunction(twoJ, twoN, twoM,  -gamma, -beta, -alpha )) );
    REQUIRE( std::conj(yap::DFunction(twoJ, twoM, twoN, alpha, beta, gamma)) == double(yap::pow_negative_one((twoN - twoM) / 2)) * yap::DFunction(twoJ, -twoM, -twoN,  alpha, beta, gamma ) );
}

TEST_CASE( "WignerD" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    // choose arbitrary angles to test with
    double alpha = 0.4 * yap::pi<double>();
    double beta = 0.6 * yap::pi<double>();
    double gamma = 0.5 * yap::pi<double>();

    SECTION( "Caching" ) {
        // cache J = 0, 1/2, 1 (out of order)
        REQUIRE_NOTHROW( yap::dMatrix::cache(0) );
        REQUIRE_NOTHROW( yap::dMatrix::cache(2) );
        REQUIRE_NOTHROW( yap::dMatrix::cache(1) );
    }

    SECTION( "J = 0") {
        SECTION ("d matrix") {
            // check invalid args
            REQUIRE_THROWS_AS( yap::dFunction(0, 1, 0, beta), yap::exceptions::Exception );
            REQUIRE_THROWS_AS( yap::dFunction(0, 0, 1, beta), yap::exceptions::Exception );

            // check val when M or N exceeds J
            REQUIRE( yap::dFunction(0, 2, 0, beta) == 0 );
            REQUIRE( yap::dFunction(0, 0, 2, beta) == 0 );

            // check val
            REQUIRE( yap::dFunction(0, 0, 0, beta) == 1 );
        }

        SECTION ("D matrix") {
            // check symmetries
            checkDSymmetries(0, 0, 0, alpha, beta, gamma);
            checkDSymmetries(0, 2, 0, alpha, beta, gamma);
            checkDSymmetries(0, 0, 2, alpha, beta, gamma);
        }
    }

    SECTION( "J = 1/2") {
        SECTION( "d matrix") {
            // check invalid args
            REQUIRE_THROWS_AS( yap::dFunction(1, 2, 0, beta), yap::exceptions::Exception );
            REQUIRE_THROWS_AS( yap::dFunction(1, 0, 2, beta), yap::exceptions::Exception );

            // check val when M or N exceeds J
            REQUIRE( yap::dFunction(1, 3, 0, beta) == 0 );
            REQUIRE( yap::dFunction(1, 0, 3, beta) == 0 );

            // check vals
            REQUIRE( yap::dFunction(1, +1, +1, beta) == Approx(+cos(beta / 2)) );
            REQUIRE( yap::dFunction(1, +1, -1, beta) == Approx(-sin(beta / 2)) );
            REQUIRE( yap::dFunction(1, -1, +1, beta) == Approx(+sin(beta / 2)) );
            REQUIRE( yap::dFunction(1, -1, -1, beta) == Approx(+cos(beta / 2)) );

            // check symmetries
            REQUIRE( yap::dFunction(1, +1, +1, beta) - yap::dFunction(1, -1, -1, beta) == 0);
            REQUIRE( yap::dFunction(1, -1, +1, beta) + yap::dFunction(1, +1, -1, beta) == 0);
            REQUIRE( yap::dFunction(1, +1, +1, beta) - yap::dFunction(1, +1, +1, -beta) == 0);
            REQUIRE( yap::dFunction(1, +1, -1, beta) - yap::dFunction(1, -1, +1, -beta) == 0);
        }

        SECTION ("D matrix") {
            // check symmetries
            checkDSymmetries(1, +1, +1, alpha, beta, gamma);
            checkDSymmetries(1, +1, -1, alpha, beta, gamma);
            checkDSymmetries(1, -1, +1, alpha, beta, gamma);
            checkDSymmetries(1, -1, -1, alpha, beta, gamma);
        }
    }

    SECTION("J = 1") {
        SECTION( "d matrix") {
            // check invalid args
            REQUIRE_THROWS_AS( yap::dFunction(2, 1, 0, beta), yap::exceptions::Exception );
            REQUIRE_THROWS_AS( yap::dFunction(2, 0, 3, beta), yap::exceptions::Exception );

            // check val when M or N exceeds J
            REQUIRE( yap::dFunction(2, 4, 0, beta) == 0 );
            REQUIRE( yap::dFunction(2, 0, 4, beta) == 0 );

            // check vals
            REQUIRE( Approx(yap::dFunction(2, +2, +2, beta)) == (1 + cos(beta)) / 2 );
            REQUIRE( Approx(yap::dFunction(2, +2,  0, beta)) == -sin(beta) / sqrt(2) );
            REQUIRE( Approx(yap::dFunction(2, +2, -2, beta)) == (1 - cos(beta)) / 2 );
            REQUIRE( Approx(yap::dFunction(2,  0, +2, beta)) == +sin(beta) / sqrt(2) );
            REQUIRE( Approx(yap::dFunction(2,  0,  0, beta)) == cos(beta) );
            REQUIRE( Approx(yap::dFunction(2,  0, -2, beta)) == -sin(beta) / sqrt(2) );
            REQUIRE( Approx(yap::dFunction(2, -2, +2, beta)) == (1 - cos(beta)) / 2 );
            REQUIRE( Approx(yap::dFunction(2, -2,  0, beta)) == +sin(beta) / sqrt(2) );
            REQUIRE( Approx(yap::dFunction(2, -2, -2, beta)) == (1 + cos(beta)) / 2 );

            // check symmetries
            REQUIRE( yap::dFunction(2, +2,  0, beta) + yap::dFunction(2, -2,  0, beta) == 0);
            REQUIRE( yap::dFunction(2, +2,  0, beta) + yap::dFunction(2,  0, +2, beta) == 0);
            REQUIRE( yap::dFunction(2, +2, -2, beta) - yap::dFunction(2, -2, +2, -beta) == 0);
            REQUIRE( yap::dFunction(2,  0, +2, beta) - yap::dFunction(2, +2,  0, -beta) == 0);
        }

        SECTION ("D matrix") {
            // check symmetries
            checkDSymmetries(2, +2, +2, alpha, beta, gamma);
            checkDSymmetries(2, +2,  0, alpha, beta, gamma);
            checkDSymmetries(2, +2, -2, alpha, beta, gamma);
            checkDSymmetries(2,  0, +2, alpha, beta, gamma);
            checkDSymmetries(2,  0,  0, alpha, beta, gamma);
            checkDSymmetries(2,  0, -2, alpha, beta, gamma);
            checkDSymmetries(2, -2, +2, alpha, beta, gamma);
            checkDSymmetries(2, -2,  0, alpha, beta, gamma);
            checkDSymmetries(2, -2, -2, alpha, beta, gamma);
        }
    }
}
