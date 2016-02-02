#include <catch.hpp>

#include <Constants.h>
#include <Exceptions.h>
#include <logging.h>
#include <WignerD.h>

#include <cmath>

TEST_CASE( "WignerD" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    // choose an arbitrary angle to test with
    double beta = 0.6 * yap::PI;

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
    }
}
