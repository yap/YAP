#include <catch.hpp>

#include <ClebschGordan.h>
#include <Exceptions.h>
#include <logging.h>

TEST_CASE( "ClebschGordan" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    SECTION( "Error throwing" ) {

        // nonzero coefficient
        REQUIRE_THROWS_AS( yap::ClebschGordan::nonzeroCoefficient(1, 1, 1, 1, 1, 0), yap::exceptions::AngularMomentumNotConserved);
        REQUIRE_THROWS_AS( yap::ClebschGordan::nonzeroCoefficient(0, 1, 1, 1, 1, 0), yap::exceptions::InconsistentSpinProjection);

        // coefficient
        REQUIRE_THROWS_AS( yap::ClebschGordan::coefficient(1, 1, 1, 1, 1, 0), yap::exceptions::AngularMomentumNotConserved);
        REQUIRE_THROWS_AS( yap::ClebschGordan::coefficient(0, 1, 1, 1, 1, 0), yap::exceptions::InconsistentSpinProjection);

        // nonzero coupling
        REQUIRE_THROWS_AS( yap::ClebschGordan::nonzeroCoupling(1, 1, 1, 1, 1, 0, 1), yap::exceptions::AngularMomentumNotConserved);
        REQUIRE_THROWS_AS( yap::ClebschGordan::nonzeroCoupling(1, 1, 1, 0, 1, 0, 0), yap::exceptions::InconsistentSpinProjection);
    }

    SECTION ( "Coefficients" ) {

        SECTION ( "Symmetries" ) {
            unsigned maxJ = 2;
            for (unsigned two_j1 = 0; two_j1 <= 2 * maxJ; ++two_j1)
                for (unsigned two_j2 = 0; two_j2 <= 2 * maxJ; ++two_j2)
                    for (int two_m1 = -two_j1; two_m1 <= (int)two_j1; two_m1 += 2)
                        for (int two_m2 = -two_j2; two_m2 <= (int)two_j2; two_m2 += 2)
                            for (unsigned two_J = std::max(abs((int)two_j1 - (int)two_j2), abs(two_m1 + two_m2)); two_J <= two_j1 + two_j2; two_J += 2) {
                                // negate m1 and m2
                                REQUIRE(
                                    yap::ClebschGordan::coefficient(two_j1, two_m1, two_j2, two_m2, two_J)
                                    - yap::pow_negative_one((two_J - two_j1 - two_j2) / 2) * yap::ClebschGordan::coefficient(two_j1, -two_m1, two_j2, -two_m2, two_J)
                                    == Approx(0)
                                );
                                // 1 <--> 2
                                REQUIRE(
                                    yap::ClebschGordan::coefficient(two_j1, two_m1, two_j2, two_m2, two_J)
                                    - yap::pow_negative_one((two_J - two_j1 - two_j2) / 2) * yap::ClebschGordan::coefficient(two_j2, two_m2, two_j1, two_m1, two_J)
                                    == Approx(0)
                                );
                            }
        }

        SECTION ( "j2 = 0" ) {
            for (unsigned two_j1 = 0; two_j1 <= 10; ++two_j1)
                for (int two_m1 = -two_j1; two_m1 <= (int)two_j1; two_m1 += 2) {
                    REQUIRE( yap::ClebschGordan::coefficient(two_j1, two_m1, 0, 0, two_j1, two_m1) == Approx( 1 ) );
                    REQUIRE( yap::ClebschGordan::coefficient(two_j1, two_m1, 0, 0, two_j1 + 2, two_m1) == Approx( 0 ) );
                    if (two_m1 < (int)two_j1)
                        REQUIRE( yap::ClebschGordan::coefficient(two_j1, two_m1, 0, 0, two_j1, two_m1 + 2) == Approx( 0 ) );
                }
        }

        SECTION ( "j1 = 1/2, j2 = 1/2" ) {
            REQUIRE( yap::ClebschGordan::coefficient(1, +1, 1, +1, 2, +2) == Approx( 1 ) );
            REQUIRE( yap::ClebschGordan::coefficient(1, -1, 1, -1, 2, -2) == Approx( 1 ) );
            REQUIRE( yap::ClebschGordan::coefficient(1, +1, 1, -1, 2, 0) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(1, +1, 1, -1, 0, 0) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(1, -1, 1, +1, 2, 0) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(1, -1, 1, +1, 0, 0) == Approx(-sqrt(1. / 2)) );
        }

        SECTION ( "j1 = 1, j2 = 1" ) {
            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2, +2, 4, 4) == Approx( 1 ) );

            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2,  0, 4, 2) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(2,  0, 2, +2, 4, 2) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2,  0, 2, 2) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(2,  0, 2, +2, 2, 2) == Approx(-sqrt(1. / 2)) );

            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2, -2, 4, 0) == Approx(sqrt(1. / 6)) );
            REQUIRE( yap::ClebschGordan::coefficient(2,  0, 2,  0, 4, 0) == Approx(sqrt(2. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(2, -2, 2, +2, 4, 0) == Approx(sqrt(1. / 6)) );
            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2, -2, 2, 0) == Approx(sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(2,  0, 2,  0, 2, 0) == Approx( 0 ) );
            REQUIRE( yap::ClebschGordan::coefficient(2, -2, 2, +2, 2, 0) == Approx(-sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(2, +2, 2, -2, 0, 0) == Approx( sqrt(1. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(2,  0, 2,  0, 0, 0) == Approx(-sqrt(1. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(2, -2, 2, +2, 0, 0) == Approx( sqrt(1. / 3)) );
        }

        SECTION ( "j1 = 2, j2 = 1" ) {
            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2, +2, 6, 6) == Approx(1) );

            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2,  0, 6, 4) == Approx( sqrt(1. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2, +2, 6, 4) == Approx( sqrt(2. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2,  0, 4, 4) == Approx( sqrt(2. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2, +2, 4, 4) == Approx(-sqrt(1. / 3)) );

            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2, -2, 6, 2) == Approx( sqrt(1. / 15)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2,  0, 6, 2) == Approx( sqrt(8. / 15)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2, +2, 6, 2) == Approx( sqrt(2. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2, -2, 4, 2) == Approx( sqrt(1. / 3)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2,  0, 4, 2) == Approx( sqrt(1. / 6)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2, +2, 4, 2) == Approx(-sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +4, 2, -2, 2, 2) == Approx( sqrt(3. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2,  0, 2, 2) == Approx(-sqrt(3. / 10)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2, +2, 2, 2) == Approx( sqrt(1. / 10)) );

            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2, -2, 6, 0) == Approx( sqrt(1. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2,  0, 6, 0) == Approx( sqrt(3. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, -2, 2, +2, 6, 0) == Approx( sqrt(1. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2, -2, 4, 0) == Approx( sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2,  0, 4, 0) == Approx( 0 ) );
            REQUIRE( yap::ClebschGordan::coefficient(4, -2, 2, +2, 4, 0) == Approx(-sqrt(1. / 2)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, +2, 2, -2, 2, 0) == Approx( sqrt(3. / 10)) );
            REQUIRE( yap::ClebschGordan::coefficient(4,  0, 2,  0, 2, 0) == Approx(-sqrt(2. / 5)) );
            REQUIRE( yap::ClebschGordan::coefficient(4, -2, 2, +2, 2, 0) == Approx( sqrt(3. / 10)) );
        }

    }
}
